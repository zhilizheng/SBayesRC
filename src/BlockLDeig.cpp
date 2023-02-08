/* SBayesRC
 *
 * This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
 * Develped by Zhili Zheng <zhilizheng@outlook.com>, 2021
 */

#include "commR.hpp"
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include "BlockLDeig.hpp"

using std::to_string;
using std::endl;

void BlockLDeig::readLD(string mldmDir, double cutThresh, const VectorXf &bhat, string outPrefix){
    string mldinfo = mldmDir + "/ldm.info";
    std::ifstream info(mldinfo.c_str());
    if(!info){
        Rcout << "Can't read " << mldinfo << ". Please check the file contains the LD information." << std::endl;
        throw("error");
    }
    vector<string> header = {"block","chr", "idxBlockStart","idxBlockEnd","idxStarts","idxEnds","preBlock","postBlock"};
    string line;
    std::getline(info, line);
    std::istringstream line_buf(line);
    std::istream_iterator<string> begin(line_buf), end;
    vector<string> line_elements(begin, end);
    if(line_elements.size() == 8){
        if(line_elements != header){
            Rcout << "Error: invalid header in " << mldinfo << std::endl;
            throw("error");
        }
    }else{
        Rcout << "Error: invalid file (incorrect format): " << mldinfo << std::endl;
        throw("error");
    }

    int line_number = 1;
    while(std::getline(info, line)){
        std::istringstream line_buf(line);
        std::istream_iterator<string> begin(line_buf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8){
            Rcout << "Invalid line " << line_number + 1 << " in " << mldinfo << std::endl;
            throw("error");
        }
        /*
        if(stoi(line_elements[0]) != line_number){
            Rcout << "block can only be in sequence, error in line: " << line_number + 1<< std::endl; 
            throw("error");
        }
        */

        idxBlocks.push_back(stoi(line_elements[0]));

        int idxstart = stoi(line_elements[2]) - 1;
        int idxend = stoi(line_elements[3]) - 1;
        startPos.push_back(idxstart);
        endPos.push_back(idxend);
        ms.push_back(idxend - idxstart + 1);

        idxStarts.push_back(stoi(line_elements[4]) - 1);
        idxEnds.push_back(stoi(line_elements[5]) - 1);
        idxBlkPres.push_back(stoi(line_elements[6]) - 1);
        idxBlkAfters.push_back(stoi(line_elements[7]) - 1);

        line_number++;
    }

    nBlocks = idxStarts.size();
    m = endPos[nBlocks-1] + 1;


    //init size
    qs = VectorXf::Zero(nBlocks);
    w.resize(nBlocks);
    Q.resize(nBlocks);

    Rcout << " Number of variants from LD information: " << m << endl;
    Rcout << "Start reading LD information, and cutting the variance to " << cutThresh << "..." << endl;

    vector<string> exts = {"var0.9999", "var0.9995", "var0.999", "var0.995", "var0.99"};

    string ext = "";
    for(int i = 0; i < exts.size(); i++){
        string testldm = mldmDir + "/eig_block1_" + exts[i] + ".bin3";
        if(FILE *fp = fopen(testldm.c_str(), "rb")){
            fclose(fp);
            ext = "_" + exts[i] + ".bin3";
        }
    }

    if(ext == ""){
        Rcout << "Error: can't find eigen binary file." << std::endl;
        throw("Error");
    }


    vector<int> cur_ms(nBlocks);
    vector<int> cur_ks(nBlocks);
    vector<int> set_ks(nBlocks);
    vector<float> w_sum(nBlocks);
    vector<float> Q_sum(nBlocks);
    #pragma omp parallel for schedule(dynamic)
    for(int idx = 0; idx < nBlocks; idx++){
        int curBlock = idxBlocks[idx];
        string ldm = mldmDir + "/eig_block" + to_string(curBlock) + ext;
        FILE *fp = fopen(ldm.c_str(), "rb");
        if(!fp){
            Rcout << "Error to read file eig_block" << curBlock << ext << std::endl;
            exit(1);
        }
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumLambda = 0;
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            Rcout << "Read " << ldm << " error (m)" << endl;
            throw("read file error");
        }
        if(cur_m != ms[idx]){
            Rcout << ldm << ": inconsistent marker number to marker information" << std::endl;
            Rcout << ms[idx] << " marker: " << cur_m << std::endl;
            throw("error");
        }

        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            Rcout << "Read " << ldm << " error (k)" << endl;
            throw("read file error");
        }

        if(fread(&sumLambda, sizeof(float), 1, fp) != 1){
            Rcout << "Read " << ldm << " error sumLambda" << endl;
            throw("read file error");
        }

        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            Rcout << "Read " << ldm << " error (lambda)" << endl;
            throw("read file error");
        }

        // cut thresh further
        float sums = 0;
        float varThresh = cutThresh * sumLambda;
        int32_t set_k = cur_k;
        for(int j = 0; j < cur_k; j++){
            sums += lambda[j];
            if(sums >= varThresh){
                set_k = j+1;
                break;
            }
        }

        if(set_k > cur_k){
            set_k = cur_k;
        }


        qs[idx] = set_k;
        MatrixXf U(cur_m, set_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)set_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            Rcout << "Read " << ldm << " error (U)" << endl;
            throw("read file error");
        }
        fclose(fp);

        VectorXf sqrtLambda = lambda.head(set_k).array().sqrt();
        VectorXf curw = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (U.transpose() * bhat.segment(startPos[idx], ms[idx]));
        MatrixXf curQ = sqrtLambda.asDiagonal() * U.transpose();
        w[idx] = curw;
        Q[idx] = curQ;

        cur_ms[idx] = cur_m;
        cur_ks[idx] = cur_k;
        set_ks[idx] = set_k;
        w_sum[idx] = curw.sum();
        Q_sum[idx] = curQ.sum();

        /*
        if((idx+1) % 100 == 0){
            double t100 = (timer.now() - tic)/1e9;
            readTime += t100;
            Rcout << " read 100 block up to " << idx << ", time: " << t100 << "..." << std::endl;
            tic = timer.now();
        }
        */
    }

    //output information
    if(!outPrefix.empty()){
        string outfile = outPrefix + ".rlog"; 
        std::ofstream fout(outfile.c_str());
        if(!fout){
            Rcout << "Error: can't write to " << outfile << std::endl;
            exit(101);
        }
        fout << "idx\tBlock\tStartPos\tEndPos\tLDsize\tLDSsize\tInUse\tw_sum\tQ_sum" << std::endl;

        for(int idx = 0; idx < nBlocks; idx++){
            fout << (idx + 1) << "\t" << idxBlocks[idx] << "\t" << startPos[idx] << "\t" << endPos[idx] << "\t" << cur_ms[idx] << "\t" << cur_ks[idx] << "\t" << set_ks[idx] << "\t" << w_sum[idx] << "\t" << Q_sum[idx] << std::endl;
        }

    }
    Rcout << "Finish reading " << nBlocks << " LD blocks" << std::endl;

    /*
       VectorXf weightQ(m);
       weightQ.fill(0);
       VectorXf weightQn(m);
       weightQn.fill(0);
       for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
       int baseStart = startPos[idxBlk];
       int baseEnd = endPos[idxBlk];
       const MatrixXf &curQ = Q[idxBlk];
       for(int idx = baseStart; idx <= baseEnd; idx++){
       VectorXf curQi = curQ.col(idx - baseStart);
       weightQ(idx) += curQi.dot(curQi);
       weightQn(idx) += 1;
       }
       }

       weightQ = weightQ.array() / weightQn.array();

       VectorXf weightQ2(m);
       for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
       int baseStart = startPos[idxBlk];
       int baseEnd = endPos[idxBlk];
       const MatrixXf &curQ = Q[idxBlk];

       int startM = idxStarts[idxBlk];
       int endM = idxEnds[idxBlk];

       for(int idx = startM; idx <= endM; idx++){
       VectorXf curQi = curQ.col(idx - baseStart);
       weightQ2(idx) = curQi.dot(curQi);
       }
       }
       */
}

int BlockLDeig::getIdxStart(int blk){
    return idxStarts[blk];
}

int BlockLDeig::getIdxEnd(int blk){
    return idxEnds[blk];
}

int BlockLDeig::getStartPos(int blk){
    return startPos[blk];
}

float BlockLDeig::getq(int blk){
    return qs[blk];
}

int BlockLDeig::getNBlocks(){
    return nBlocks;
}

int BlockLDeig::getNMarker(){
    return m;
}

Ref<const MatrixXf> BlockLDeig::getQ(int blk) const{
    return Q[blk];
}

Ref<VectorXf> BlockLDeig::getW(int blk){
    return w[blk];
}


