/* Timer
 *
 * This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
 * Develped by Zhili Zheng <zhilizheng@outlook.com>, 2021
 */

#include "Timer.h"
#include "commR.h"

std::map<string, std::chrono::time_point<std::chrono::steady_clock>> Timer::time_map;

void Timer::start(std::string key){
    time_map[key] = std::chrono::steady_clock::now();
}

float Timer::elapse(string key){
    auto end = std::chrono::steady_clock::now();
    float secs = 0.0;
    if(time_map.find(key) != time_map.end()){
        auto duration = end - time_map[key];
        secs = std::chrono::duration_cast<std::chrono::duration<float>>(duration).count();
    }else{
        Rcout << " Error: can't find key of time start" << std::endl;
        throw("error");
    }
    return secs;
}


