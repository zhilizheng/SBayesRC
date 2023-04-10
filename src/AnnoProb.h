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

#ifndef SBRC_ANNO_PROB_H
#define SBRC_ANNO_PROB_H

#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

using namespace std;

using Eigen::Map;
using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::MatrixXf;
using Eigen::MatrixXd;
using Eigen::Matrix2f;

class AnnoProb {
   private:
    MatrixXf annoMat;
    VectorXf annoSD;
    VectorXf XPXiSD;
    vector<Matrix2f> XPXqiSD;
    MatrixXf annoDists;
    MatrixXf vg_annot_comp;
    VectorXf vg_tot_annot;
    VectorXf vg_enrich_qt;
    VectorXf vg_enrich_bin;
    bool bOutDetail;

    vector<bool> bAnnoBinary;

    map<string, ofstream> outs;

    VectorXf invFrac;
    MatrixXf alpha;
    MatrixXf annoCondProb;
    MatrixXf annoJointProb;

    MatrixXf snpP;
    VectorXf sigmaSq;
    VectorXf ssq;
    int numDist;
    int numComp;
    int numAnno;
    int numSnps;
    void initP_Pi_annoAlpha(const VectorXf &pi, MatrixXf &snpPi);
    void annoCondProb_probit();
    void annoSigmaSS_sample();
    void annoEffect_sample_Gibbs(const MatrixXf &z);
    void compute_AnnoJointProb();
    void computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi);
    void writeHeader(const vector<string> &annoStrs);

   public:
    AnnoProb(string fileAnnot, int numAnno,  const VectorXf &Pi, MatrixXf &snpPi, bool bOutDetail, double initSS=1);
    void samplePi(const MatrixXf &z, MatrixXf &snpPi);
    void computeEnrichBin(const MatrixXf &vg_snp_comp, float vp);
    void computeEnrichQt(const VectorXf &vg_snp, float vg);

    void computeDist(const MatrixXf &z, const VectorXf &numSnpMix);

    void computeProb();

    void open(string prefix, const vector<string> &annoStrs);
    void close();
    void output();
    VectorXd getSigmaSq();
};

#endif //SBRC_ANNO_PROB_H

