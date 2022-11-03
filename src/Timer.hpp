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

#ifndef TIMER_H
#define TIMER_H
#include <chrono>
#include <map>
#include <string>

using std::string;

class Timer{
public:
    void start(string maker);
    float elapse(string marker);

private:
    static std::map<string, std::chrono::time_point<std::chrono::steady_clock>> time_map;

};
#endif //TIMMER_H

