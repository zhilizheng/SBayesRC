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

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>

#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <Rcpp/Benchmark/Timer.h>


using namespace Rcpp;
using namespace std;
using namespace arma;

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::MatrixXf;

//' cut the LD with new thresh
//' @param ldm string, input ldm
//' @param outfile string, output file
//' @cutThresh double, cut threshold
//' @export
// [[Rcpp::export]]
bool cutLD(std::string ldm, std::string outfile, double cutThresh = 1){
    Timer timer;
    nanotime_t tic = timer.now();

    FILE *fp = fopen(ldm.c_str(), "rb");
    if(!fp){
        Rcout << "Error to read LD file" << std::endl;
        throw("read file error!");
    }

    int32_t cur_m = 0;
    int32_t cur_k = 0;
    float sumLambda = 0;

    if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
        Rcout << "Read " << ldm << " error (m)" << endl;
        throw("read m error");
    }

    if (fread(&cur_k, sizeof(int32_t), 1, fp) != 1) {
        Rcout << "Read " << ldm << " error (k)" << endl;
        throw("read file error");
    }

    if (fread(&sumLambda, sizeof(float), 1, fp) != 1) {
        Rcout << "Read " << ldm << " error sumLambda" << endl;
        throw("read file error");
    }

    VectorXf lambda(cur_k);
    if (fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k) {
        Rcout << "Read " << ldm << " error (lambda)" << endl;
        throw("read file error");
    }

    // cut thresh further
    float sums = 0;
    float varThresh = cutThresh * sumLambda;
    int32_t set_k = cur_k;
    for (int j = 0; j < cur_k; j++) {
        sums += lambda[j];
        if (sums >= varThresh) {
            set_k = j + 1;
            break;
        }
    }

    if (set_k > cur_k) {
        set_k = cur_k;
    }

    MatrixXf U(cur_m, set_k);
    uint64_t nElements = (uint64_t)cur_m * (uint64_t)set_k;
    if (fread(U.data(), sizeof(float), nElements, fp) != nElements) {
        Rcout << "Read " << ldm << " error (U)" << endl;
        throw("read file error");
    }
    fclose(fp);

    Rcout << "Prepare and reading time: " << ((timer.now() - tic) / 1e9 ) << endl;
    tic = timer.now();

    VectorXf curLambda = lambda.head(set_k);

    
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

    if(fwrite(curLambda.data(), sizeof(float), set_k, oFile) != set_k){
        Rcout << "Write " << outfile << " error lambda" << endl;
        throw("write error");
    }

    if (fwrite(U.data(), sizeof(float), nElements, oFile) != nElements) {
        Rcout << "Write " << ldm << " error (U)" << endl;
        throw("write error");
    }
    fclose(oFile);

    Rcout << "Finished" << endl;
    return true;

}
