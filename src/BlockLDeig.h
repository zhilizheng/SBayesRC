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

#ifndef BLOCK_LD_EIG_H
#define BLOCK_LD_EIG_H
#include <vector>
#include <Eigen/Eigen>
#include <string>

using std::vector;
using Eigen::VectorXf;
using Eigen::MatrixXf;
using Eigen::Ref;
using std::string;

class BlockLDeig{
public:
    void readLD(string mldm, double cutTresh, const VectorXf &bhat, string outPrefix);

    int getIdxStart(int blk);
    int getIdxEnd(int blk);
    int getStartPos(int blk);
    float getq(int blk);
    int getNBlocks();
    int getNMarker();

    Ref<const MatrixXf> getQ(int blk) const;
    Ref<VectorXf> getW(int blk);

private:
    int nBlocks; // number of blocks
    int m; // number of marker

    vector<int> idxBlocks;
    vector<int> idxStarts, idxEnds;
    vector<int> idxBlkPres, idxBlkAfters;

    vector<int> startPos;
    vector<int> endPos;
    vector<int> ms;

    VectorXf qs; 
    vector<VectorXf> w;
    vector<MatrixXf> Q;
};

int getLDPrefix(string mldm, double &cutThresh, string &tempstr); 
bool read1LD(string tempstr, int idxBlock, int type, float cutThresh, MatrixXf &U, VectorXf &lambda, float &sumLambda);

#endif //BLOCK_LD_EIG_H

