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

#ifndef SBRC_DIST_H
#define SBRC_DIST_H

#include <Eigen/Eigen>

using Eigen::VectorXf;

namespace InvChiSq{
    float sample(const float df, const float scale);
}

namespace Bernoulli{
    int sample(const VectorXf &p);
}


namespace Gamma{
    float sample(const float shape, const float scale);
}


namespace Dirichlet{
    VectorXf sample(const int n, const VectorXf &irx);
}


namespace Normal{
    double quantile_01(double value);
    double cdf_01(double value);
    double sample(double mean, double variance);
}


namespace TruncatedNormal {
    double sample_tail_01_rejection(double a);
    double sample_upper_truncated(double mean, double sd, double b);
    double sample_lower_truncated(double mean, double sd, double a);
}

#endif //SBRC_DIST_H
