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

#include "commR.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include "BlockLDeig.h"

using std::to_string;
using std::endl;

void BlockLDeig::readLD(string mldmDir, double cutThresh, const VectorXf &bhat, string outPrefix){
    string mldinfo = mldmDir + separator() + "ldm.info";
    std::ifstream info(mldinfo.c_str());
    if(!info){
        Rcout << "Can't read " << mldinfo << ". Please check the file contains the LD information." << std::endl;
        throw("error");
    }
    vector<string> header1 = {"block","chr", "idxBlockStart","idxBlockEnd","idxStarts","idxEnds","preBlock","postBlock"};
    vector<string> header2 = {"Block", "Chrom", "StartSnpIdx", "StartSnpID", "EndSnpIdx", "EndSnpID", "NumSnps"};
    string line;
    std::getline(info, line);
    std::istringstream line_buf(line);
    std::istream_iterator<string> begin(line_buf), end;
    vector<string> line_elements(begin, end);
    bool bValidInfo = false;
    if(line_elements.size() == 8){
        if(line_elements == header1){
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
            bValidInfo = true;
        }
    }

    if(line_elements.size() == 7){
        if(line_elements == header2){
            int line_number = 1;
            while(std::getline(info, line)){
                std::istringstream line_buf(line);
                std::istream_iterator<string> begin(line_buf), end;
                vector<string> line_elements(begin, end);
                if(line_elements.size() != 7){
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

                int idxstart = stoi(line_elements[2]);
                int idxend = stoi(line_elements[4]);
                startPos.push_back(idxstart);
                endPos.push_back(idxend);
                int curM = stoi(line_elements[6]);
                if(curM != idxend - idxstart + 1){
                    Rcout << "Inconsistent number of markers in line " <<  line_number + 1 << " in " << mldinfo << std::endl;
                    throw("Error");
                }
                ms.push_back(curM);

                idxStarts.push_back(idxstart);
                idxEnds.push_back(idxend);
                idxBlkPres.push_back(-1);
                idxBlkAfters.push_back(-1);

                line_number++;
            }
            bValidInfo = true;
        }
    }
 
    if(!bValidInfo){
        Rcout << "Error: invalid file (incorrect format): " << mldinfo << std::endl;
        throw("error");
    }

    nBlocks = idxStarts.size();
    m = endPos[nBlocks-1] + 1;

    //init size
    qs = VectorXf::Zero(nBlocks);
    w.resize(nBlocks);
    Q.resize(nBlocks);

    Rcout << " Number of variants from LD information: " << m << endl;
    Rcout << "Start reading LD information, and cutting the variance to " << cutThresh << "..." << endl;

    string tempstr;
    double fileThresh;
    int type = getLDPrefix(mldmDir, fileThresh, tempstr);
    if(type < 1){
        Rcout << "can't find the valid LD file in folder " << mldmDir << std::endl;
        throw("Error");
    }else{
        Rcout << "Found LD in " << tempstr << ", type: " << type << std::endl;
    }

    if(fileThresh < cutThresh){
        Rcout << "can't set cutting threshold smaller than threshold in file, " << "Threshold in file: " << fileThresh << ", treshold set: " << cutThresh << std::endl;
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
        VectorXf lambda(1);
        MatrixXf U(1,1);
        float sumLambda;

        bool status = read1LD(tempstr, curBlock, type, cutThresh, U, lambda, sumLambda);
        if(!status){
            Rcout << "can't read file with string template " << tempstr << std::endl;
            throw("Error");
        }

        if(U.rows() != ms[idx]){
            Rcout << tempstr << " block " << curBlock << " inconsistent marker number to marker information" << std::endl;
            Rcout << ms[idx] << " marker: " << U.rows() << std::endl;
            throw("error");
        }

        VectorXf sqrtLambda = lambda.array().sqrt();
        VectorXf curw = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (U.transpose() * bhat.segment(startPos[idx], ms[idx]));
        MatrixXf curQ = sqrtLambda.asDiagonal() * U.transpose();
        w[idx] = curw;
        Q[idx] = curQ;

        cur_ms[idx] = U.rows();
        set_ks[idx] = U.cols();
        qs[idx] = U.cols();
        w_sum[idx] = curw.sum();
        Q_sum[idx] = curQ.sum();
    }

    //output information
    if(!outPrefix.empty()){
        string outfile = outPrefix + ".rlog"; 
        std::ofstream fout(outfile.c_str());
        if(!fout){
            Rcout << "Error: can't write to " << outfile << std::endl;
            throw("Error");
        }
        fout << "idx\tBlock\tStartPos\tEndPos\tLDsize\tInUse\tw_sum\tQ_sum" << std::endl;

        for(int idx = 0; idx < nBlocks; idx++){
            fout << (idx + 1) << "\t" << idxBlocks[idx] << "\t" << startPos[idx] << "\t" << endPos[idx] << "\t" << cur_ms[idx] << "\t" << set_ks[idx] << "\t" << w_sum[idx] << "\t" << Q_sum[idx] << std::endl;
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

// get the LDprefix return the version
// 0, don't find 
// 1, old version, with bin3
// 2. new version with eigen.bin
int getLDPrefix(string mldm, double &cutThresh, string &tempstr){
    vector<string> exts = {"0.9995", "0.999", "0.995", "0.99"};

    string curExt = "";
    int retCode = 0;
    for(int i = 0; i < exts.size(); i++){
        string testldm = mldm + separator() + "eig_block1_var" + exts[i] + ".bin3";
        if(FILE *fp = fopen(testldm.c_str(), "rb")){
            fclose(fp);
            curExt =  exts[i];
            retCode = 1;
            tempstr = mldm + separator() + "eig_block{BLOCK}_var" + exts[i] + ".bin3";
            cutThresh = std::stod(curExt); 
        }
    }

    if(curExt == ""){
        string testldm = mldm + separator() + "block1.eigen.bin";
        if(FILE *fp = fopen(testldm.c_str(), "rb")){
            float fCutThresh = 0;
            fseek(fp, 12, SEEK_SET);
            if(fread(&fCutThresh, sizeof(float), 1, fp) != 1){
                Rcout << "Read " << testldm << " error thresh" << endl;
            }else{
                retCode = 2;
                tempstr = mldm + separator() + "block{BLOCK}.eigen.bin";
                cutThresh = fCutThresh;
            }
 
            fclose(fp);
       }
    }

    return retCode;
}

#ifndef _STAND_ALONE_
//' find the eigen ld file
//' @param ldm string, input ldm
//' @return list, type, LD threshold and each LD template
// @export
// [[Rcpp::export]]
List getLDPrefix(std::string mldm){
    string tempstr;
    int retCode;
    double cutThresh;
    retCode = getLDPrefix(mldm, cutThresh, tempstr);
    return(List::create(_["type"] = retCode, 
                _["thresh"] = cutThresh,
                _["template"] = tempstr
                ));
}

// [[Rcpp::export]]
SEXP getPseudoRand(std::string tempstr, int type, int m, Rcpp::NumericVector blocks, Rcpp::NumericVector startPos, double thresh, Eigen::Map<Eigen::VectorXd> rand){
    VectorXf bt_add(m);
    VectorXf randf = rand.cast<float>();
    int n = blocks.size();
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n; i++){
        int curBlock = blocks[i];
        MatrixXf U(1,1);
        VectorXf lambda(1);
        float sumLambda;
        bool status = read1LD(tempstr, curBlock, type, thresh, U, lambda, sumLambda);
        if(!status){
            Rcout << "Read block " << curBlock << " error" << std::endl;
            throw("Error");
        }
        bt_add.segment(startPos[i], U.rows()) = U * (lambda.array().sqrt().matrix().asDiagonal() * randf.segment(startPos[i], U.rows()));
    }

    return Rcpp::wrap(bt_add.cast<double>());
}
#endif 

#include <boost/algorithm/string.hpp>
bool read1LD(string tempstr, int idxBlock, int type, float cutThresh, MatrixXf &U, VectorXf &lambda, float &sumLambda){
    boost::replace_all(tempstr, "{BLOCK}", std::to_string(idxBlock));
    string ldm = tempstr;
    FILE *fp = fopen(ldm.c_str(), "rb");
    if(!fp){
        Rcout << "Read " << ldm << " error (m)" << endl;
        throw("Error");
    }
    int32_t cur_m = 0;
    int32_t cur_k = 0;

    if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
        Rcout << "Read " << ldm << " error (m)" << endl;
        throw("read file error");
    }
    if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
        Rcout << "Read " << ldm << " error (k)" << endl;
        throw("read file error");
    }

    if(fread(&sumLambda, sizeof(float), 1, fp) != 1){
        Rcout << "Read " << ldm << " error sumLambda" << endl;
        throw("read file error");
    }

    if(type == 2){
        float fileThresh = 0;
        if(fread(&fileThresh, sizeof(float), 1, fp) != 1){
            Rcout << "Read " << ldm << " error thresh" << endl;
            throw("read file error");
        }
    }

    lambda.resize(cur_k);
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


    lambda.conservativeResize(set_k);
    U.resize(cur_m, set_k);
    uint64_t nElements = (uint64_t)cur_m * (uint64_t)set_k;
    if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
        Rcout << "Read " << ldm << " error (U)" << endl;
        throw("read file error");
    }
    fclose(fp);
    return true;
}
