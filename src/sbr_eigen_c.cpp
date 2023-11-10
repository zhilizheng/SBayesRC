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

#include "SBayesRC.h"

#include <Eigen/Eigen>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace std;

using Eigen::Map;
using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::MatrixXf;
using Eigen::MatrixXd;


// [[Rcpp::export]]
List sbayesr_eigen_joint_annot(int niter, int burn, Eigen::Map<Eigen::VectorXd> bhat, int numAnno, Rcpp::StringVector annoStrs, std::string mldmDir, double vary, Eigen::Map<Eigen::VectorXd> blkN, Eigen::Map<Eigen::VectorXd> cgamma, Eigen::Map<Eigen::VectorXd> startPi, Rcpp::IntegerVector rmSNPIndices, double starth2=0.01, double cutThresh=1, bool bOrigin = false, std::string outPrefix="", std::string samVe = "fixVe", double resam_thresh=1.1, bool bOutDetail=false, int outFreq=10, double initAnnoSS=1.0, bool bOutBeta=false){

    string curSamVe = samVe;
    VectorXf fbhat = bhat.cast<float>();
    VectorXf n = blkN.cast<float>();
    VectorXf fgamma = cgamma.cast<float>();
    VectorXf pi = startPi.cast<float>();

    vector<string> annoStrings(annoStrs.size());
    for(int i = 0; i < annoStrs.size(); i++){
        annoStrings[i] = annoStrs[i];
    }

    vector<int> v_rmSNPIndices(rmSNPIndices.size());
    for(int i = 0; i < rmSNPIndices.size(); i++){
        v_rmSNPIndices[i] = rmSNPIndices[i];
    }

    SBayesRC sbr(niter, burn, fbhat, numAnno, annoStrings, mldmDir, vary, n, fgamma, pi, v_rmSNPIndices, starth2, cutThresh, bOrigin, outPrefix, samVe, resam_thresh, bOutDetail, outFreq, initAnnoSS);
    sbr.setOutBeta(bOutBeta);
    sbr.mcmc();

    sbr.outRmIndex(outPrefix);

    // return 
    VectorXd mean_par_vec = sbr.get_mean_par_vec();
    NumericVector mean_par(mean_par_vec.data(), mean_par_vec.data() + mean_par_vec.size());
    mean_par.attr("names") = CharacterVector::create("hsq", "nnz", "sigmaSq", "ssq", "vare", "varg");

    VectorXf betaMean_vec = sbr.get_betaMean_vec();
    NumericVector betaMean(betaMean_vec.data(), betaMean_vec.data()+betaMean_vec.size());

    VectorXf hsq_mcmc = sbr.get_hsq_mcmc();
    VectorXf hsq2_mcmc = sbr.get_hsq2_mcmc();
    NumericVector hsq_mcmc_r(hsq_mcmc.data(), hsq_mcmc.data() + hsq_mcmc.size());
    NumericVector hsq2_mcmc_r(hsq2_mcmc.data(), hsq2_mcmc.data() + hsq2_mcmc.size());

    VectorXf beta = sbr.get_beta();
    NumericVector beta_last_r(beta.data(), beta.data() + beta.size());

    MatrixXf pi_mcmc = sbr.get_pi_mcmc();
    NumericMatrix pi_mcmc_r(pi_mcmc.rows(), pi_mcmc.cols(), pi_mcmc.data());

    MatrixXf vg_comp_mcmc = sbr.get_vg_comp_mcmc();
    NumericMatrix vg_comp_mcmc_r(vg_comp_mcmc.rows(), vg_comp_mcmc.cols(), vg_comp_mcmc.data());

    MatrixXf n_comp_mcmc = sbr.get_n_comp_mcmc();
    NumericMatrix n_comp_mcmc_r(n_comp_mcmc.rows(), n_comp_mcmc.cols(), n_comp_mcmc.data());

    MatrixXf pip = sbr.get_pip();
    NumericMatrix pip_r(pip.rows(), pip.cols(), pip.data());

    //NumericVector weightQ_r(weightQ.data(), weightQ.data() + weightQ.size());
    Eigen::VectorXd sigmaAnno = sbr.get_anno_ss();
    NumericVector sigmaAnno_r(sigmaAnno.data(), sigmaAnno.data() + sigmaAnno.size());
 
    return(List::create(_["par"]=mean_par, 
                _["betaMean"] = betaMean, 
                _["betaLast"] = beta_last_r,
                _["hsq_hist"] = hsq_mcmc_r, 
                _["ssq_hist"] = hsq2_mcmc_r, 
                _["pi_hist"] = pi_mcmc_r,
                _["vare_hist"] = sbr.get_vare_infos(),
                _["ssq_block_hist"] = sbr.get_ssq_infos(),
                _["hsq_block_hist"] = sbr.get_hsq_infos(),
                _["pip"] = pip_r,
                _["vg_comp_hist"] = vg_comp_mcmc_r,
                _["n_comp_hist"] = n_comp_mcmc_r,
                _["sigma_anno"] = sigmaAnno_r
                )); 
}
