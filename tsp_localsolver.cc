

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <ostream>
#include <fstream>





#include "ortools/base/commandlineflags.h"
// #include "absl/flags/flag.h"
#include "localsolver.h"
// #include "tsptw_data_dt.h"
#include "./localsolver_result.pb.h"
#include "./localsolver_vrp.pb.h"

DEFINE_string(instance_file, "", "instance file name or data");
DEFINE_string(solution_file, "", "solution file name ");

DEFINE_int64(time_limit_in_ms, 1000, "Time limit in ms, no option means no limit.");
DEFINE_int64(no_solution_improvement_limit, -1, "Iterations whitout improvement");
DEFINE_int64(minimum_duration, -1, "Initial time whitout improvement in ms");
DEFINE_int64(init_duration, -1, "Maximum duration to find a first solution");
DEFINE_int64(time_out_multiplier, 2, "Multiplier for the nexts time out");
DEFINE_int64(vehicle_limit, 0, "Define the maximum number of vehicle");
DEFINE_int64(solver_parameter, -1, "Force a particular behavior");
DEFINE_bool(only_first_solution, false, "Compute only the first solution");
DEFINE_bool(verification_only, false, "Only verify the suplied initial solution");
DEFINE_bool(balance, false, "Route balancing");
DEFINE_bool(nearby, false, "Short segment priority");
#ifdef DEBUG
DEFINE_bool(debug, true, "debug display");
#else
DEFINE_bool(debug, false, "debug display");
#endif
DEFINE_bool(intermediate_solutions, false, "display intermediate solutions");
DEFINE_string(routing_search_parameters,
              "", /* An example of how we can override the default settings */
              // "first_solution_strategy:ALL_UNPERFORMED"
              // "local_search_operators {"
              // "  use_path_lns:BOOL_TRUE"
              // "  use_inactive_lns:BOOL_TRUE"
              // "}",
              "Text proto RoutingSearchParameters (possibly partial) that will "
              "override the DefaultRoutingSearchParameters()");

using namespace localsolver;
using namespace std;

int main(int argc, char **argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  gflags::ParseCommandLineFlags(&argc, &argv, true);
   

  localsolver_vrp::Problem problem;

  const std::string filename = absl::GetFlag(FLAGS_instance_file);
  std::fstream input(filename, ios::in | ios::binary);
  if (!problem.ParseFromIstream(&input)) {
    cout << "Failed to parse pbf." << endl;
  }
  // LocalSolver localsolver;
  // LSModel model = localsolver.getModel(); 
  


  // for (auto vehicle : problem.vehicles()){

  //   cout << "vehicle : " << vehicle.id() << endl;
  // }

  

  for (auto vehicle : problem.vehicles())
  {
    vehicle.PrintDebugString();
  }
  


  //=================================================================

  //=================================================================


  gflags::ShutDownCommandLineFlags();

  return 0;
}

