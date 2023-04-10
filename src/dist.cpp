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

#include "dist.h"

#include "commR.h"
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <random>


using namespace std;


namespace InvChiSq{
    typedef boost::mt19937 random_engine;
    typedef boost::gamma_distribution<> gamma_distribution;
    typedef boost::variate_generator<random_engine&, gamma_distribution> gamma_generator;
    static thread_local random_engine engine;
    float sample(const float df, const float scale){
        gamma_generator sgamma(engine, gamma_distribution(0.5f*df, 1));
        return scale/(2.0f*sgamma());
    }
}

namespace Bernoulli{
    typedef boost::uniform_01<> uniform_01;
    typedef boost::mt19937 random_engine;
    typedef boost::variate_generator<random_engine&, uniform_01> uniform01_generator;
    //static random_engine engine;
    //static uniform01_generator ranf(engine, uniform_01());

    static thread_local std::mt19937_64 rng;
    static thread_local std::uniform_real_distribution<float> unif(0, 1);

    int sample(const VectorXf &p){
        float cum = 0;
        //float rnd = ranf();
        float rnd = unif(rng);
        long size = p.size();
        unsigned ret = 0;
        for (unsigned i=0; i<size; ++i) {
            if (!isnan(p[i])) cum += p[i];
            if (rnd < cum) {
                ret = i;
                break;
            }
        }
        return ret;
    }
}

namespace Gamma{
    typedef boost::mt19937 random_engine;
    typedef boost::gamma_distribution<> gamma_distribution;
    typedef boost::variate_generator<random_engine&, gamma_distribution> gamma_generator;
    static thread_local random_engine engine;
    float sample(const float shape, const float scale){
        gamma_generator sgamma(engine, gamma_distribution(shape, scale));
        return sgamma();
    }
}

namespace Dirichlet{
    VectorXf sample(const int n, const VectorXf &irx){
        VectorXf ps(n);
        double sx = 0.0;
        for(int i = 0; i < n; i++){
            ps[i] = Gamma::sample(irx(i), 1.0);
            sx += ps[i];
        }
        ps = ps / sx;
        return ps;
    }
}

namespace Normal{
    typedef boost::mt19937 random_engine;
    static thread_local random_engine engine;
    typedef boost::uniform_01<> uniform_01;
    typedef boost::normal_distribution<> normal_distribution;
    typedef boost::variate_generator<random_engine&, uniform_01> uniform01_generator;
    typedef boost::variate_generator<random_engine&, normal_distribution> normal_generator;
    static thread_local uniform01_generator ranf(engine, uniform_01());
    static thread_local normal_generator snorm(engine, normal_distribution(0,1));  // standard normal

    boost::math::normal_distribution <> d = boost::math::normal_distribution <> (0 ,1);
    double quantile_01(double value){
        // qnorm(value)
        return quantile(d, value);
    }

    double cdf_01(double value){
        // pnorm(value)
        return cdf(d, value);
    }

    double sample(double mean, double variance){
        return mean + snorm()*sqrt(variance);
    }
}

namespace TruncatedNormal {

    double sample_tail_01_rejection(double a) {  // a < x < inf
        double b = 0.5 * a * a;
        double w, v;
        do {
            double u = Normal::ranf();
            w = b - log(u);
            v = Normal::ranf();
        } while (v > w / b);
        return sqrt(2.0 * w);
    }

    double sample_lower_truncated(double mean, double sd, double a) {  // a < x < inf
        //    int seed = rand();
        //    float x = truncated_normal_a_sample (mean, sd, a, seed);
        //    return x;

        if (a - mean > 5 * sd) {
            return mean + sd * sample_tail_01_rejection((a - mean) / sd);
        } else {
            double alpha = (a - mean) / sd;
            double alpha_cdf = Normal::cdf_01(alpha);
            double u = Normal::ranf();
            double x = alpha_cdf + (1.0 - alpha_cdf) * u;  // x ~ Uniform(alpha_cdf, 1)
            //if (x <= 0 || x >= 1) Rcout << "alpha " << alpha << " alpha_cdf " << alpha_cdf << " u " << u << " x " << x << std::endl;
            return mean + sd * Normal::quantile_01(x);
        }
    }

    double sample_upper_truncated(double mean, double sd, double b) {  // -inf < x < b
        //    int seed = rand();
        //    float x = truncated_normal_b_sample (mean, sd, b, seed);
        //    return x;
        if (mean - b > 5.0 * sd) {
            return mean - sd * sample_tail_01_rejection((mean - b) / sd);
        } else {
            double beta = (b - mean) / sd;
            double beta_cdf = Normal::cdf_01(beta);
            double u;
            do {
                u = Normal::ranf();
            } while (!u);
            double x = beta_cdf * u;  // x ~ Uniform(0, beta_cdf);
            //if (x <= 0 || x >= 1) Rcout << "beta " << beta << " beta_cdf " << beta_cdf << " u " << u << " x " << x << std::endl;
            return mean + sd * Normal::quantile_01(x);
        }
    }

}  // namespace TruncatedNormal

/*
NumericVector sample_upper(Map<VectorXd> means, Map<VectorXd> sds, Map<VectorXd> bs){
    NumericVector ret(means.size());
    for(int i = 0; i < means.size(); i++){
        ret[i] = TruncatedNormal::sample_upper_truncated(means[i], sds[i], bs[i]);
    }
    return(ret);
}

NumericVector sample_lower(Map<VectorXd> means, Map<VectorXd> sds, Map<VectorXd> bs){
    NumericVector ret(means.size());
    for(int i = 0; i < means.size(); i++){
        ret[i] = TruncatedNormal::sample_lower_truncated(means[i], sds[i], bs[i]);
    }
    return(ret);
}


*/
