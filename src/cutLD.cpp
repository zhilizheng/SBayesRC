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

#ifdef _OPENMP
  #include <omp.h>
#endif

//#define ARMA_64BIT_WORD 1
//#include <RcppArmadillo.h>

#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "Timer.h"
#include "commR.h"
#include "BlockLDeig.h"

using namespace Rcpp;
using namespace std;
//using namespace arma;

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::MatrixXf;

// [[Rcpp::export]]
bool cutLDc(std::string tempstr, int type, Rcpp::NumericVector blocks, std::string outDir, double cutThresh = 0.995){
    Timer timer;
    timer.start("cutLD");

    int n = blocks.size();
    float thresh = cutThresh;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n; i++){
        int curBlock = blocks[i];
        MatrixXf U(1,1);
        VectorXf lambda(1);
        float sumLambda;
        bool status = read1LD(tempstr, curBlock, type, cutThresh, U, lambda, sumLambda);
        if(!status){
            Rcout << "Read block " << curBlock << " error" << std::endl;
            throw("Error");
        }
        int32_t cur_m = U.rows();
        int32_t set_k = U.cols();
        string outfile = outDir + separator() + "block" + std::to_string(curBlock) + ".eigen.bin";
        FILE * oFile = fopen(outfile.c_str(), "wb");
        if(fwrite(&cur_m, sizeof(int32_t), 1, oFile) != 1){
            Rcout << "Write " << outfile << " error m" << endl;
            throw("write error");
        }

        if(fwrite(&set_k, sizeof(int32_t), 1, oFile) != 1){
            Rcout << "Write " << outfile << " error k" << endl;
            throw("write error");
        }

        if(fwrite(&sumLambda, sizeof(float), 1, oFile) != 1){
            Rcout << "Write " << outfile << " error sumLambda" << endl;
            throw("write error");
        }

        if(fwrite(&thresh, sizeof(float), 1, oFile) != 1){
            Rcout << "Write " << outfile << " error thresh" << endl;
            throw("write error");
        }

        if(fwrite(lambda.data(), sizeof(float), set_k, oFile) != set_k){
            Rcout << "Write " << outfile << " error lambda" << endl;
            throw("write error");
        }

        uint64_t nElements = (uint64_t) cur_m * (uint64_t) set_k;
        if (fwrite(U.data(), sizeof(float), nElements, oFile) != nElements) {
            Rcout << "Write " << outfile << " error (U)" << endl;
            throw("write error");
        }
        fclose(oFile);
    }

    Rcout << "Finished" << endl;
    return true;
}
