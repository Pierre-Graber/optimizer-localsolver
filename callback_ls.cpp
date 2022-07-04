#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <lsarray.h>
#include <lscollection.h>
#include <lsexpression.h>
#include <lssolution.h>
#include <ostream>
#include <string>
#include <typeinfo>
#include <vector>

#include "./localsolver_result.pb.h"
#include "./localsolver_vrp.pb.h"
#include "localsolver.h"
#include "ortools/base/commandlineflags.h"

DEFINE_string(solution_file, "", "solution file name ");
DEFINE_string(instance_file, "", "instance file name or data");

DEFINE_int64(time_limit_in_ms, 1000, "Time limit in ms, no option means no limit.");
DEFINE_int64(no_solution_improvement_limit, -1, "Iterations whitout improvement");
DEFINE_int64(minimum_duration, -1, "Initial time whitout improvement in ms");
DEFINE_int64(init_duration, -1, "Maximum duration to find a first solution");
DEFINE_int64(time_out_multiplier, 2, "Multiplier for the next time out");
DEFINE_int64(vehicle_limit, 0, "Define the maximum number of vehicle");
DEFINE_int64(solver_parameter, -1, "Force a particular behavior");
DEFINE_bool(only_first_solution, false, "Compute only the first solution");
DEFINE_bool(verification_only, false, "Only verify the suplied initial solution");
DEFINE_bool(balance, false, "Route balancing");
DEFINE_bool(nearby, false, "Short segment priority");

using namespace localsolver;
using namespace std;
class MyCallback : public LSCallback {
private:
  int lastBestRunningTime;
  lsint lastBestValue;

public:
  MyCallback() {
    lastBestRunningTime = 0;
    lastBestValue = 0;
  }

  void callback(LocalSolver& ls, LSCallbackType type) {
    (void)type;
    LSStatistics stats = ls.getStatistics();
    LSExpression obj = ls.getModel().getObjective(0);

    if (obj.getValue() > lastBestValue) {
      lastBestRunningTime = stats.getRunningTime();
      lastBestValue = obj.getValue();
    }

    if (stats.getRunningTime() > absl::GetFlag(FLAGS_minimum_duration) &&
        stats.getRunningTime() >
            std::max(absl::GetFlag(FLAGS_time_out_multiplier) * lastBestRunningTime,
                     absl::GetFlag(FLAGS_time_limit_in_ms))) {
      cout << ">>>>>>>  Resolution stopped" << endl;
      ls.stop();
    } else {
      cout << ">>>>>> Objective improved by " << (obj.getValue() - lastBestValue) << endl;
    }
  }
};
