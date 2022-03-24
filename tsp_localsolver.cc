

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

using namespace localsolver;
using namespace google::protobuf;
using namespace std;

class localsolver_VRP {
public:
  const localsolver_vrp::Problem problem;

  LocalSolver localsolver;
  LSModel model;
  vector<LSExpression> timeMatrices;

  LSExpression serviceMatrixIndicies;
  LSExpression serviceTime;

  LSExpression vehicleStartIndicies;
  LSExpression vehicleEndIndicies;

  void MatrixBuilder(vector<LSExpression>& Matrices, const RepeatedField<float>& matrix) {
    LSExpression Matrix(model.array());
    int matrix_size = sqrt(matrix.size());
    vector<int> zeros(matrix_size + 1, 0);
    for (int i = 0; i < matrix_size; i++) {
      vector<int> matrix_row;
      matrix_row.reserve(matrix_size);
      for (int j = 0; j < matrix_size; j++) {
        matrix_row.push_back(ceil(matrix.at(i * matrix_size + j)));
      }
      matrix_row.push_back(0);
      Matrix.addOperands(model.array(matrix_row.begin(), matrix_row.end()));
    }
    Matrix.addOperands(model.array(zeros.begin(), zeros.end()));
    Matrices.push_back(Matrix);
  }

  // Constructor
  localsolver_VRP(localsolver_vrp::Problem& pb)
      : problem(pb)
      , localsolver()
      , model(localsolver.getModel())
      , timeMatrices(0, model.array()) {
    for (const auto& matrix : problem.matrices()) {
      MatrixBuilder(timeMatrices, matrix.time());
    }

    vector<int> vehicleStartIndicies_vec;
    vector<int> vehicleEndIndicies_vect;
    for (auto vehicle : problem.vehicles()) {
      int matrix_size =
          std::max(sqrt(problem.matrices(vehicle.matrix_index()).time_size()),
                   sqrt(problem.matrices(vehicle.matrix_index()).distance_size()));
      if (vehicle.start_index() == -1) {
        vehicleStartIndicies_vec.push_back(matrix_size);
      } else {
        vehicleStartIndicies_vec.push_back(vehicle.start_index());
      }
      if (vehicle.end_index() == -1) {
        vehicleEndIndicies_vect.push_back(matrix_size);
      } else {
        vehicleEndIndicies_vect.push_back(vehicle.end_index());
      }
    }
    vehicleStartIndicies =
        model.array(vehicleStartIndicies_vec.begin(), vehicleStartIndicies_vec.end());
    vehicleEndIndicies =
        model.array(vehicleEndIndicies_vect.begin(), vehicleEndIndicies_vect.end());

    vector<int> serviceTime_vec;
    vector<int> serviceMatrixIndicies_vec;
    for (auto service : problem.services()) {
      serviceMatrixIndicies_vec.push_back(service.matrix_index());
      serviceTime_vec.push_back(service.duration());
    }
    serviceMatrixIndicies =
        model.array(serviceMatrixIndicies_vec.begin(), serviceMatrixIndicies_vec.end());
    serviceTime = model.array(serviceTime_vec.begin(), serviceTime_vec.end());
  }

  void createModelAndSolve() {
    // Decision Variables :

    // A tour for each vehicle
    vector<LSExpression> serviceSequences;
    serviceSequences.reserve(problem.vehicles_size());
    for (auto vehicle : problem.vehicles()) {
      serviceSequences.push_back(model.listVar(problem.services_size()));
    }

    // Are the vehicles actually used?
    vector<LSExpression> vehiclesUsed;

    // Number of vehicles used
    LSExpression nbVehiclesUsed;
    vehiclesUsed.resize(problem.vehicles_size());

    vector<LSExpression> routeDuration(problem.vehicles_size());
    vector<LSExpression> endTime(problem.vehicles_size());

    // all services must be satisfied by the vehicles
    LSExpression partition =
        model.partition(serviceSequences.begin(), serviceSequences.end());
    model.constraint(partition);

    vector<LSExpression> endTimes(problem.vehicles_size());

    // Constraints for each vehicle
    int k = 0;
    for (const auto& vehicle : problem.vehicles()) {
      LSExpression sequence = serviceSequences[k];
      sequence.setName("sequence_" + to_string(k));
      LSExpression c = model.count(sequence);

      vehiclesUsed[k] = c > 0;
      vehiclesUsed[k].setName("vehicleUsed_" + to_string(k));

      LSExpression timeSelector = model.createLambdaFunction([&](LSExpression i,
                                                                 LSExpression prev) {
        return model.iif(
            c == 0, 0,
            model.iif(i == 0,
                      model.at(timeMatrices[vehicle.matrix_index()],
                               vehicleStartIndicies[k],
                               serviceMatrixIndicies[sequence[0]]) +
                          serviceTime[sequence[i]],
                      prev + model.iif(i < c,

                                       model.at(timeMatrices[vehicle.matrix_index()],
                                                serviceMatrixIndicies[sequence[i - 1]],
                                                serviceMatrixIndicies[sequence[i]]) +
                                           serviceTime[sequence[i]],
                                       model.at(timeMatrices[vehicle.matrix_index()],
                                                serviceMatrixIndicies[sequence[c - 1]],
                                                vehicleEndIndicies[k]))));
      });
      routeDuration[k] = model.array(model.range(0, c + 1), timeSelector);
      endTimes[k] = routeDuration[k][c];
      ;
      // model.constraint(endTimes[k] <= vehicle.duration());
      k++;
    }

    // Total vehicles used :
    nbVehiclesUsed = model.sum(vehiclesUsed.begin(), vehiclesUsed.end());

    // Total distance traveled :
    LSExpression totalDuration;
    totalDuration = model.sum(endTimes.begin(), endTimes.end());
    // totalDuration.setName("totalDuration");

    // model.minimize(nbVehiclesUsed);
    model.minimize(totalDuration);

    model.close();

    localsolver.getParam().setTimeLimit(10);

    // With ls a LocalSolver object
    auto iis = localsolver.computeInconsistency();
    std::cout << iis.toString() << std::endl;

    localsolver.solve();

    for (int v = 0; v < problem.vehicles_size(); v++) {
      if (vehiclesUsed[v].getValue() != 1)
        continue;
      cout << serviceSequences[v].getCollectionValue().toString() << endl;
      cout << "total duration = " << totalDuration.getValue() << endl;
    }
  }
};

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  localsolver_vrp::Problem problem;

  const std::string filename = absl::GetFlag(FLAGS_instance_file);
  std::fstream input(filename, ios::in | ios::binary);
  if (!problem.ParseFromIstream(&input)) {
    cout << "Failed to parse protobuf." << endl;
  }

  for (auto vehicle : problem.vehicles()) {
    vehicle.PrintDebugString();
  }

  for (auto matrix : problem.matrices()) {
    for (int i = 0; i < sqrt(matrix.time_size()); i++) {
      for (int j = 0; j < sqrt(matrix.time_size()); j++) {
        cout << matrix.time(i * sqrt(matrix.time_size()) + j);
      }
      cout << endl;
    }
  }

  for (auto service : problem.services())
    cout << "duration of service " << service.id() << " : " << service.duration() << endl;

  localsolver_VRP model(problem);
  model.createModelAndSolve();

  gflags::ShutDownCommandLineFlags();

  return 0;
}
