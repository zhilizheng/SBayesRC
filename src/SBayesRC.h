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

#ifndef SBAYESRC_HPP
#define SBAYESRC_HPP
#include <Eigen/Eigen>
#include <vector>
#include <set>

#include "AnnoProb.h"
#include "BlockLDeig.h"

using Eigen::VectorXf;
using Eigen::MatrixXd;
using Eigen::MatrixXf;

class SBayesRC{

public:
    SBayesRC(int niter, int burn, VectorXf fbhat, int numAnno, std::vector<string> &annoStrs, std::string mldmDir, double vary, VectorXf n, VectorXf fgamma, VectorXf pi, const std::vector<int> &rmSNPIndices, double starth2=0.01, double cutThresh=1, bool bOrigin = false, std::string outPrefix="", std::string samVe = "fixVe", double resam_thresh=1.1, bool bOutDetail=false, int outFreq=10, double initAnnoSS=1.0);
    void mcmc();
    VectorXd get_mean_par_vec();
    VectorXf get_betaMean_vec();
    VectorXf get_betaMean2_vec();
    VectorXf get_betaMean3_vec();
    VectorXf get_hsq_mcmc();
    VectorXf get_hsq2_mcmc();
    VectorXf get_beta();
    MatrixXf get_pi_mcmc();
    MatrixXf get_n_comp_mcmc();
    MatrixXf get_vg_comp_mcmc();
    MatrixXf get_pip();
    MatrixXf get_pip2();
    MatrixXf get_pip3();
    MatrixXd get_vare_infos();
    MatrixXd get_hsq_infos();
    MatrixXd get_ssq_infos();
    VectorXd get_anno_ss();
    void outRmIndex(std::string outPrefix);
    void setOutFreq();
    void setOutBeta(bool bOut);

private:
    int ndist;  // number of mixture distribution components
    int m; // number of marker
    int nBlocks; // number of blocks
    VectorXf vare; // residual
    float varg; // Genetic variance
    float vary; // phenotypic variance
    VectorXf fgamma; // gamma
    bool estimateSigmaSq = true;
    bool bAnnot; // annot or not
    MatrixXf snpPi; // SNPs pi
    VectorXf pi;
    bool bOutDetail;
    string outPrefix;
    string curSamVe;
    double resam_thresh;
    VectorXf n;
    AnnoProb * anno = NULL; // annoation
    VectorXd annoSS;
    BlockLDeig blockLDeig; // block LD eigen


    VectorXf b;
    float betaThresh; 
    std::vector<std::set<int>> delSNPs;
    VectorXf beta; // joint beta
    VectorXf betasum; // sum of beta
    VectorXf betasum_all; // sum of beta all iter
    VectorXf betasum2; // sum of beta2
    VectorXf betasum3; // sum of beta3
    bool bOutBeta = false;
                      //
    int niter;
    int burn;

    bool bOrigin;

    float sigmaSq; 

    // mcmc history
    VectorXf hsq_mcmc;
    VectorXf hsq2_mcmc;
    MatrixXf pi_mcmc;
    MatrixXf n_comp_mcmc;
    MatrixXf vg_comp_mcmc;

    // out put history
    int outFreq = 10;
    MatrixXd iter_infos;
    MatrixXd vare_infos;
    MatrixXd ssq_infos;
    MatrixXd hsq_infos;
    MatrixXf pip_count;
    MatrixXf pip_count2;
    MatrixXf pip_count3;
};

#endif  //SBAYESRC_HPP
