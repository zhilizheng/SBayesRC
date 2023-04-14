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

#ifndef SBRC_ALTR_H
#define SBRC_ALTR_H

#ifdef _OPENMP
  #include <omp.h>
#endif

// use the Rcout in cpp code
#ifdef _STAND_ALONE_
#include <iostream>
#define Rcout std::cout
#include <cstdio>
#define Rprintf printf

#else
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#endif

inline char separator(){
#if defined _WIN32 || defined __CYGWIN__
    return '\\';
#else
    return '/';
#endif
}

#endif //SBRC_ALTR_H
