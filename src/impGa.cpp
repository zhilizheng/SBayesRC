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
//using Eigen::LLT;

// [[Rcpp::export]]
NumericVector impGa(std::string ldm, Eigen::Map<Eigen::VectorXd> z, Eigen::Map<Eigen::VectorXi> typedIndex, int m, double cutThresh = 1, double diag_mod = 0.1){
    Timer timer;
    nanotime_t tic = timer.now();
    VectorXf fz = z.cast<float>();
    int n_tt = typedIndex.size();
    if(n_tt != fz.size()){
        Rcout << "the index and the bhat shall be in the same size" << endl;
        throw("error");
    }
    int n_imp = m - n_tt;
    if(n_imp < 0){
        Rcout << "m is smaller than typed SNPs, impossible" << endl;
        throw("error");
    }
    
    if(n_imp == 0){
        Rcout << "Don't need to impute" << endl;
        NumericVector retImp(n_imp);
        return(retImp);
    }

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

    if(cur_m != m){
        Rcout << "m is inconsistent in file,  provided: " << m << ", in LD: " << cur_m << endl;
        throw("error");
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

    Rcout << "LD reconstruct with m = " << cur_m << ", k = " << set_k << ", full_k = " << cur_k << std::endl;
    Rcout << "Numer of SNP typed: " << n_tt << ", Number of SNP to impute: " << n_imp << std::endl;
    MatrixXf LD = U * curLambda.asDiagonal() * U.transpose();
    //MatrixXf LD = U * U.transpose();
    U.resize(0, 0);
    LD.diagonal().array() += (float)diag_mod;


    Rcout << "Reconstruct LD time: " << ((timer.now() - tic) / 1e9 ) << endl;
    tic = timer.now();

    MatrixXf LDtt(n_tt, n_tt);
    for(int i = 0; i < n_tt; i++){
        int ori_i = typedIndex[i] - 1;
        for(int j = 0; j < n_tt; j++){
            int ori_j = typedIndex[j] - 1;
            LDtt(j, i) = LD(ori_j, ori_i);
        }
    }


    Rcout << "Get LDtt time: " << ((timer.now() - tic) / 1e9 ) << endl;
    tic = timer.now();

    vector<int> full_index(m);
    iota(full_index.begin(), full_index.end(), 1);
    vector<int> remain_index(n_imp);
    set_difference(full_index.begin(), full_index.end(), 
                   typedIndex.data(), typedIndex.data() + typedIndex.size(),
                   remain_index.begin());

    MatrixXf LDit(n_imp, n_tt);
    for(int i = 0; i < n_tt; i++){
        int ori_i = typedIndex[i] - 1;
        for(int j = 0; j < n_imp; j++){
            int ori_j = remain_index[j] - 1;
            LDit(j, i) = LD(ori_j, ori_i);
        }
    }

    Rcout << "Get LDit time: " << ((timer.now() - tic) / 1e9 ) << endl;
    tic = timer.now();
    LD.resize(0, 0);

    //LLT<MatrixXf> llt;
    //llt.compute(LDtt);

    //Rcout << "LLT time: " << ((timer.now() - tic) / 1e9 ) << endl;
    //tic = timer.now();

    //VectorXf LDi_Z = llt.solve(fz);
    fmat LDtt_mat(LDtt.data(), n_tt, n_tt, false, true);
    fvec fz_vec(fz.data(), n_tt, false, true);

    fvec LDi_Z = solve(LDtt_mat, fz_vec,  solve_opts::fast + solve_opts::likely_sympd);

    Rcout << "LDi_Z time: " << ((timer.now() - tic) / 1e9 ) << endl;
    tic = timer.now();

    VectorXf LDi_Z_eig = Eigen::Map<VectorXf>(LDi_Z.memptr(), n_tt);

    VectorXf imp_z = LDit * LDi_Z_eig;

    Rcout << "imp_z time: " << ((timer.now() - tic) / 1e9 ) << endl;

    NumericVector imp_z_r(imp_z.data(), imp_z.data() + imp_z.size());
    return(imp_z_r);
}
