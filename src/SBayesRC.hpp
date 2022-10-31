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

#include "AnnoProb.hpp"
#include "BlockLDeig.hpp"

using Eigen::VectorXf;
using Eigen::MatrixXd;
using Eigen::MatrixXf;

class SBayesRC{

public:
    SBayesRC(int niter, int burn, VectorXf fbhat, int numAnno, vector<string> &annoStrs, std::string mldmDir, double vary, VectorXf n, VectorXf fgamma, VectorXf pi, double starth2=0.01, double cutThresh=1, bool bOrigin = false, std::string outPrefix="", std::string samVe = "fixVe", double resam_thresh=1.1, bool bOutDetail=false);
    void mcmc();
    VectorXd get_mean_par_vec();
    VectorXf get_betaMean_vec();
    VectorXf get_hsq_mcmc();
    VectorXf get_hsq2_mcmc();
    VectorXf get_beta();
    MatrixXf get_pi_mcmc();
    MatrixXf get_n_comp_mcmc();
    VectorXf get_vg_comp_mcmc();
    MatrixXf get_pip();
    MatrixXd get_vare_infos();
    MatrixXd get_hsq_infos();
    MatrixXd get_ssq_infos();

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
    BlockLDeig blockLDeig; // block LD eigen


    VectorXf beta; // joint beta
    VectorXf betasum; // sum of beta
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
};

#endif  //SBAYESRC_HPP
