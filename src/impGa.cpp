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
//

#if !defined _WIN32 && !defined __CYGWIN__
#define EIGEN_USE_BLAS
#endif
#include "commR.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "Timer.h"
#include "BlockLDeig.h"


using namespace Rcpp;
using namespace std;
//using namespace arma;

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::MatrixXf;
using Eigen::LLT;

// [[Rcpp::export]]
NumericVector impGa(std::string tempstr, int curBlock, int type, Eigen::Map<Eigen::VectorXd> z, Eigen::Map<Eigen::VectorXi> typedIndex, int m, double cutThresh = 0.995, double diag_mod = 0.1){
    Timer timer;
    timer.start("imp");
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

    MatrixXf U(1,1);
    VectorXf curLambda(1);
    float sumLambda;

    bool status = read1LD(tempstr, curBlock, type, cutThresh, U, curLambda, sumLambda);
    if(!status){
        Rcout << "Read block " << curBlock << " error" << std::endl;
        throw("Error");
    }

    Rcout << "Prepare and reading time: " << timer.elapse("imp") << endl;
    timer.start("imp");

    Rcout << "LD reconstruct with m = " << U.rows() << ", k = " << U.cols() << std::endl;
    Rcout << "Numer of SNP typed: " << n_tt << ", Number of SNP to impute: " << n_imp << std::endl;
    MatrixXf LD = U * curLambda.asDiagonal() * U.transpose();
    //MatrixXf LD = U * U.transpose();
    U.resize(0, 0);
    LD.diagonal().array() += (float)diag_mod;


    Rcout << "Reconstruct LD time: " << timer.elapse("imp") << endl;
    timer.start("imp");

    MatrixXf LDtt(n_tt, n_tt);
    for(int i = 0; i < n_tt; i++){
        int ori_i = typedIndex[i] - 1;
        for(int j = 0; j < n_tt; j++){
            int ori_j = typedIndex[j] - 1;
            LDtt(j, i) = LD(ori_j, ori_i);
        }
    }


    Rcout << "Get LDtt time: " << timer.elapse("imp") << endl;
    timer.start("imp");

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

    Rcout << "Get LDit time: " << timer.elapse("imp") << endl;
    timer.start("imp");
    LD.resize(0, 0);

    LLT<MatrixXf> llt;
    llt.compute(LDtt);

    //Rcout << "LLT time: " << ((timer.now() - tic) / 1e9 ) << endl;
    //tic = timer.now();

    VectorXf LDi_Z_eig = llt.solve(fz);
    //fmat LDtt_mat(LDtt.data(), n_tt, n_tt, false, true);
    //fvec fz_vec(fz.data(), n_tt, false, true);

    //fvec LDi_Z = solve(LDtt_mat, fz_vec,  solve_opts::fast + solve_opts::likely_sympd);

    Rcout << "LDi_Z time: " <<  timer.elapse("imp") << endl;
    timer.start("imp");

    //VectorXf LDi_Z_eig = Eigen::Map<VectorXf>(LDi_Z.memptr(), n_tt);

    VectorXf imp_z = LDit * LDi_Z_eig;

    Rcout << "imp_z time: " <<  timer.elapse("imp") << endl;

    NumericVector imp_z_r(imp_z.data(), imp_z.data() + imp_z.size());
    return(imp_z_r);
}
