

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

#define CUSTOM_MAX_INT (int64) std::pow(2, 30)

DEFINE_string(solution_file, "", "solution file name ");
DEFINE_string(instance_file, "", "instance file name or data");

using namespace localsolver;
using namespace google::protobuf;
using namespace std;

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

class localsolver_VRP {
public:
  const localsolver_vrp::Problem problem;

  LocalSolver localsolver;
  LSModel model;
  vector<LSExpression> timeMatrices; // the matrices are reordered by service.matrix_index

  LSExpression serviceTime;
  LSExpression serviceQuantitiesMatrix;
  LSExpression serviceExclusionCost;

  LSExpression vehicleStartIndicies;
  LSExpression vehicleEndIndicies;
  LSExpression vehicleCapacitiesMatrix;

  map<string, int64> ids_map_;
  vector<map<int, LSExpression>> timesFromWarehouses;
  vector<map<int, LSExpression>> timesToWarehouses;

  LSExpression twStartsArray;
  LSExpression twEndsArray;
  LSExpression twAbsoluteEndsArray;

  LSExpression nextStart(const LSExpression& service, const LSExpression& time) {
    LSExpression timeSelector = model.createLambdaFunction([&](LSExpression k) {
      return model.iif(model.at(twAbsoluteEndsArray, service, k) >= time,
                       model.at(twStartsArray, service, k),
                       model.at(twAbsoluteEndsArray, service,
                                twAbsoluteEndsArray[service].getNbOperands()));
    });

    LSExpression earliestAvailableTime = model.min(
        model.range(0, twAbsoluteEndsArray[service].getNbOperands()), timeSelector);

    return model.max(time, earliestAvailableTime);
  }

  void MatrixBuilder(vector<LSExpression>& Matrices, const RepeatedField<float>& matrix) {
    LSExpression Matrix(model.array());
    int matrix_size = sqrt(matrix.size());

    for (const auto& serviceFrom : problem.services()) {
      vector<int> matrix_row;
      matrix_row.reserve(matrix_size);
      for (const auto& serviceTo : problem.services()) {
        matrix_row.push_back(ceil(matrix.at(serviceFrom.matrix_index() * matrix_size +
                                            serviceTo.matrix_index())));
      }
      Matrix.addOperands(model.array(matrix_row.begin(), matrix_row.end()));
    }
    Matrices.push_back(Matrix);
  }

  int64 IdIndex(const string& id) const {
    map<string, int64>::const_iterator it = ids_map_.find(id);
    if (it != ids_map_.end())
      return it->second;
    else
      return -1;
  }

  string IndexId(const int index) const {
    for (auto it = ids_map_.begin(); it != ids_map_.end(); ++it) {
      if (it->second == index) {
        return it->first;
      }
    }
    return "not found";
  }

  // Constructor
  localsolver_VRP(localsolver_vrp::Problem& pb)
      : problem(pb)
      , localsolver()
      , model(localsolver.getModel())
      , timeMatrices(0, model.array())
      , serviceQuantitiesMatrix(model.array())
      , vehicleCapacitiesMatrix(model.array())
      , twStartsArray(model.array())
      , twEndsArray(model.array())
      , twAbsoluteEndsArray(model.array()) {
    for (const auto& matrix : problem.matrices()) {
      MatrixBuilder(timeMatrices, matrix.time());
      cout << timeMatrices.rbegin()->toString() << endl;
    }

    vector<int> vehicleStartIndiciesVec;
    vector<int> vehicleEndIndiciesVec;
    vector<vector<float>> vehicleCapacitiesMatrixVec;
    int k = 0;
    timesFromWarehouses.resize(problem.matrices_size());
    timesToWarehouses.resize(problem.matrices_size());
    for (const auto& vehicle : problem.vehicles()) {
      int time_matrix_size = sqrt(problem.matrices(vehicle.matrix_index()).time_size());
      map<int, LSExpression>& vehicleStartMap =
          timesFromWarehouses[vehicle.matrix_index()];
      if (vehicleStartMap.find(vehicle.start_index()) == vehicleStartMap.end()) {
        vector<int> timesFromWarehouse;
        timesFromWarehouse.reserve(problem.services_size());
        if (vehicle.start_index() > -1) {
          for (const auto& service : problem.services()) {
            timesFromWarehouse.push_back(
                problem.matrices(vehicle.matrix_index())
                    .time(vehicle.start_index() * time_matrix_size +
                          service.matrix_index()));
          }
        } else {
          timesFromWarehouse.resize(problem.services_size());
          fill_n(timesFromWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleStartMap[vehicle.start_index()] =
            model.array(timesFromWarehouse.begin(), timesFromWarehouse.end());
      }
      map<int, LSExpression>& vehicleEndMap = timesToWarehouses[vehicle.matrix_index()];
      if (vehicleEndMap.find(vehicle.end_index()) == vehicleEndMap.end()) {
        vector<int> timesToWarehouse;
        timesToWarehouse.reserve(problem.services_size());
        if (vehicle.end_index() > -1) {
          for (const auto& service : problem.services()) {
            timesToWarehouse.push_back(
                problem.matrices(vehicle.matrix_index())
                    .time(service.matrix_index() * time_matrix_size +
                          vehicle.end_index()));
          }
        } else {
          timesToWarehouse.resize(problem.services_size());
          fill_n(timesToWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleEndMap[vehicle.end_index()] =
            model.array(timesToWarehouse.begin(), timesToWarehouse.end());
      }
      cout << vehicle.id() << " " << vehicleStartMap[vehicle.start_index()].toString()
           << endl;
      cout << vehicle.id() << " " << vehicleEndMap[vehicle.end_index()].toString()
           << endl;

      int matrix_size =
          max(sqrt(problem.matrices(vehicle.matrix_index()).time_size()),
              sqrt(problem.matrices(vehicle.matrix_index()).distance_size()));
      if (vehicle.start_index() == -1) {
        vehicleStartIndiciesVec.push_back(matrix_size);
      } else {
        vehicleStartIndiciesVec.push_back(vehicle.start_index());
      }
      if (vehicle.end_index() == -1) {
        vehicleEndIndiciesVec.push_back(matrix_size);
      } else {
        vehicleEndIndiciesVec.push_back(vehicle.end_index());
      }
      vector<float> capacityUnit;
      vehicleCapacitiesMatrixVec.push_back(capacityUnit);
      for (int unit_i = 0; unit_i < vehicle.capacities_size(); unit_i++) {
        vehicleCapacitiesMatrixVec[k].push_back(vehicle.capacities(unit_i).limit());
      }
      vehicleCapacitiesMatrix.addOperand(model.array(
          vehicleCapacitiesMatrixVec[k].begin(), vehicleCapacitiesMatrixVec[k].end()));
      k++;
    }

    vehicleStartIndicies =
        model.array(vehicleStartIndiciesVec.begin(), vehicleStartIndiciesVec.end());
    vehicleEndIndicies =
        model.array(vehicleEndIndiciesVec.begin(), vehicleEndIndiciesVec.end());

    vector<int> serviceTimeVec;
    vector<vector<float>> serviceQuantitiesVec;
    vector<float> serviceExclusionCostVec;

    int s = 0;
    for (const auto& service : problem.services()) {
      serviceTimeVec.push_back(service.duration());
      vector<float> serviceQuantityUnit;
      serviceQuantitiesVec.push_back(serviceQuantityUnit);
      for (int unit_i = 0; unit_i < service.quantities_size(); unit_i++) {
        serviceQuantitiesVec[s].push_back(service.quantities(unit_i));
      }
      serviceQuantitiesMatrix.addOperand(
          model.array(serviceQuantitiesVec[s].begin(), serviceQuantitiesVec[s].end()));
      if (service.exclusion_cost() > 0) {
        serviceExclusionCostVec.push_back(service.exclusion_cost());
      } else {
        serviceExclusionCostVec.push_back(100);
      }
      vector<int> serviceTWStarts;
      vector<int> serviceTWEnds;
      vector<int> serviceTWAbsoluteEnds;
      for (const auto& tw : service.time_windows()) {
        serviceTWStarts.push_back(tw.start());
        serviceTWEnds.push_back(tw.end());
        serviceTWAbsoluteEnds.push_back(tw.end() + tw.maximum_lateness());
      }
      serviceTWAbsoluteEnds.back();
      twStartsArray.addOperand(
          model.array(serviceTWStarts.begin(), serviceTWStarts.end()));
      twStartsArray.setName("TW starts :");
      twEndsArray.addOperand(model.array(serviceTWEnds.begin(), serviceTWEnds.end()));
      twAbsoluteEndsArray.addOperand(
          model.array(serviceTWAbsoluteEnds.begin(), serviceTWAbsoluteEnds.end()));
      ids_map_[(string)service.id()] = s;
      s++;
    }

    serviceTime = model.array(serviceTimeVec.begin(), serviceTimeVec.end());
    serviceExclusionCost =
        model.array(serviceExclusionCostVec.begin(), serviceExclusionCostVec.end());
  }

  //   void precedenceConstraints(LSModel& model) {}

  void createModelAndSolve() {
    // Decision Variables :

    // A tour for each vehicle
    vector<LSExpression> serviceSequences;
    serviceSequences.reserve(problem.vehicles_size() + 1);
    for (auto vehicle : problem.vehicles()) {
      serviceSequences.push_back(model.listVar(problem.services_size()));
    }

    // extra_super_vehicle
    serviceSequences.push_back(model.listVar(problem.services_size()));

    // Are the vehicles actually used?
    vector<LSExpression> vehiclesUsed;

    // Number of vehicles used
    LSExpression nbVehiclesUsed;
    vehiclesUsed.resize(problem.vehicles_size() + 1);

    vector<LSExpression> routeDuration(problem.vehicles_size());
    vector<LSExpression> endTime(problem.vehicles_size());
    vector<LSExpression> beginTime(problem.vehicles_size());
    vector<LSExpression> serviceStartsInTW(problem.vehicles_size());

    // all services must be satisfied by the vehicles
    LSExpression partition =
        model.partition(serviceSequences.begin(), serviceSequences.end());
    model.constraint(partition);

    // Constraints for each vehicle
    int k = 0;
    for (const auto& vehicle : problem.vehicles()) {
      LSExpression sequence = serviceSequences[k];
      sequence.setName("sequence_" + to_string(k));
      LSExpression c = model.count(sequence);

      vehiclesUsed[k] = c > 0;
      vehiclesUsed[k].setName("vehicleUsed_" + to_string(k));

      // Duration of travel of truck k
      LSExpression durationSelector = model.createLambdaFunction([&](LSExpression i) {
        return model.at(timeMatrices[vehicle.matrix_index()], sequence[i - 1],
                        sequence[i]) +
               serviceTime[sequence[i - 1]];
      });

      routeDuration[k] =
          model.sum(model.range(1, c), durationSelector) +
          model.iif(c > 0,
                    ((vehicle.has_time_window())
                         ? static_cast<lsint>(vehicle.time_window().start())
                         : 0) +
                        timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                           [sequence[0]] +
                        timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                         [sequence[c - 1]] +
                        serviceTime[sequence[c - 1]],
                    0);
      if (vehicle.duration() > 0 &&
          (!vehicle.has_time_window() ||
           vehicle.time_window().end() - vehicle.time_window().start() >
               vehicle.duration())) {
        LSExpression duration_constraint =
            (routeDuration[k] - static_cast<lsint>(vehicle.time_window().start())) <=
            static_cast<lsint>(vehicle.duration());
        duration_constraint.setName("duration_constraint_" + to_string(k));
        model.constraint(duration_constraint);
      }

      if (vehicle.has_time_window() && vehicle.time_window().end() < CUSTOM_MAX_INT) {
        model.constraint(routeDuration[k] <=
                         static_cast<lsint>(vehicle.time_window().end() +
                                            vehicle.time_window().maximum_lateness()));
      }

      // End of each visit
      LSExpression beginSelector =
          model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
            return model.min(
                model.at(twAbsoluteEndsArray, sequence[i],
                         twAbsoluteEndsArray[sequence[i]].getNbOperands() - 1),
                model.max(
                    model.at(twStartsArray, sequence[i], 0),
                    model.iif(i == 0,
                              ((vehicle.has_time_window())
                                   ? static_cast<lsint>(vehicle.time_window().start())
                                   : 0) +
                                  timesFromWarehouses[vehicle.matrix_index()]
                                                     [vehicle.start_index()][sequence[0]],
                              prev + serviceTime[sequence[i - 1]] +
                                  model.at(timeMatrices[vehicle.matrix_index()],
                                           sequence[i - 1], sequence[i]))));
          });

      beginTime[k] = model.array(model.range(0, c), beginSelector);

      LSExpression serviceTWVerifs = model.createLambdaFunction([&](LSExpression i) {
        LSExpression serviceTWVerif = model.or_(true);
        for (int tw_index = 0;
             tw_index < twAbsoluteEndsArray[sequence[i]].getNbOperands(); ++tw_index) {
          serviceTWVerif = serviceTWVerif ||
                           (model.at(beginTime[k], sequence[i]) >=
                                model.at(twStartsArray, sequence[i], tw_index) &&
                            model.at(beginTime[k], sequence[i]) <=
                                model.at(twAbsoluteEndsArray, sequence[i], tw_index));
        }

        return serviceTWVerif;
      });
      model.constraint(model.and_(model.range(0, c), serviceTWVerifs));

      // model.constraint(model.or_(model.range(0,
      // problem.services(0).time_windows_size()),
      //                            servicesTWVerifSelector));

      //   // vector<LSExpression> routeQuantitiesVec;
      //   LSExpression routeQuantities(model.array());
      for (int unit_i = 0; unit_i < vehicle.capacities_size(); unit_i++) {
        LSExpression quantityCumulator = model.createLambdaFunction([&](LSExpression i) {
          return model.at(serviceQuantitiesMatrix, sequence[i], unit_i);
        });
        LSExpression routeQuantities = model.sum(model.range(0, c), quantityCumulator);
        model.constraint(routeQuantities <= model.at(vehicleCapacitiesMatrix, k, unit_i));
      }

      for (const auto& relation : problem.relations()) {
        if (relation.type() == "shipment") {
          for (int unit = 0; unit < vehicle.capacities_size(); unit++) {
            LSExpression quantityCumulator =
                model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
                  return prev + model.at(serviceQuantitiesMatrix, sequence[i], unit);
                });
            LSExpression routeQuantityUnit =
                model.array(model.range(0, c), quantityCumulator);
            LSExpression quantityUnitChecker =
                model.createLambdaFunction([&](LSExpression i) {
                  return routeQuantityUnit[i] <= vehicle.capacities(unit).limit();
                });
            model.constraint(model.and_(model.range(0, c), quantityUnitChecker));
          }
        }
        if (relation.type() == "order" || relation.type() == "shipment") {
          for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
               link_index++) {
            model.constraint(
                model.contains(sequence, static_cast<lsint>(
                                             IdIndex(relation.linked_ids(link_index)))) ==
                model.contains(sequence, static_cast<lsint>(IdIndex(
                                             relation.linked_ids(link_index + 1)))));
            model.constraint(
                model.indexOf(sequence, static_cast<lsint>(
                                            IdIndex(relation.linked_ids(link_index)))) <=
                model.indexOf(sequence, static_cast<lsint>(IdIndex(
                                            relation.linked_ids(link_index + 1)))));
          }
        }
      }

      k++;
    }

    LSExpression unassignedServices = serviceSequences[problem.vehicles_size()];
    LSExpression c = model.count(unassignedServices);
    vehiclesUsed[problem.vehicles_size()] = c > 0;

    LSExpression exclusionCostCumulator = model.createLambdaFunction(
        [&](LSExpression i) { return serviceExclusionCost[unassignedServices[i]]; });

    LSExpression totalExclusionCost =
        model.sum(model.range(0, c), exclusionCostCumulator);
    totalExclusionCost.setName("cost");

    // Total vehicles used :
    nbVehiclesUsed = model.sum(vehiclesUsed.begin(), vehiclesUsed.end());

    // Total distance traveled :
    LSExpression totalDuration;
    totalDuration = model.sum(routeDuration.begin(), routeDuration.end());
    totalDuration.setName("totalDuration");

    // model.minimize(nbVehiclesUsed);
    model.minimize(totalExclusionCost);
    model.minimize(totalDuration);

    model.close();

    cout << model.toString() << endl;
    localsolver.getParam().setTimeLimit(5);

    // With ls a LocalSolver object
    auto iis = localsolver.computeInconsistency();
    cout << iis.toString() << endl;

    localsolver.solve();

    cout << endl;
    cout << "total duration = " << totalDuration.getValue() << endl;
    cout << "vehicle capacities matrix"
         << vehicleCapacitiesMatrix.getArrayValue().toString() << endl;

    cout << "service Quantities Matrix "
         << serviceQuantitiesMatrix.getArrayValue().toString() << endl;
    if (vehiclesUsed[problem.vehicles_size()].getValue() == 0) {
      cout << "No unassigned service" << endl;
    } else {
      cout << "unassigned service(s) : ";
      LSCollection servicesCollection =
          serviceSequences[problem.vehicles_size()].getCollectionValue();
      for (lsint i = 0; i < servicesCollection.count(); i++) {
        cout << IndexId(servicesCollection[i]) << " ";
      }
      cout << endl;
    }

    for (int v = 0; v < problem.vehicles_size(); v++) {
      if (vehiclesUsed[v].getValue() != 1)
        continue;
      cout << "assigned service(s) to " << problem.vehicles(v).id() << " : ";
      LSCollection servicesCollection = serviceSequences[v].getCollectionValue();
      for (lsint i = 0; i < servicesCollection.count(); i++) {
        cout << IndexId(servicesCollection[i]) << " ";
      }
      cout << endl;
      //   cout << "assigned service(s) to " << problem.vehicles(v).id()
      //        << " " + serviceSequences[v].getCollectionValue().toString() << endl;
    }
    for (int i = 0; i < model.getNbObjectives(); i++) {
      cout << model.getObjective(i).toString() << endl;
    }
  }
};

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  localsolver_vrp::Problem problem;

  const string filename = absl::GetFlag(FLAGS_instance_file);
  fstream input(filename, ios::in | ios::binary);
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

  for (auto service : problem.services()) {
    cout << "duration of service " << service.id() << " : " << service.duration() << endl;
  }
  localsolver_VRP model(problem);
  model.createModelAndSolve();

  gflags::ShutDownCommandLineFlags();

  return 0;
}
// // Gives the right timewindow to compare to the actual time
// LSExpression compare = model.createFunction([&](LSExpression i, LSExpression prev) {
//   return model.iif(
//       i == 0,
//       model.iif(model.at(earliestArray2, sequence[i]) == -2,
//                 model.at(earliestArray1, sequence[i]),
//                 model.iif(distanceWarehousesArray[sequence[0]] + StartDepot[k] >
//                               model.at(latestArray1, sequence[i]) -
//                                   model.at(serviceArray, sequence[i]),
//                           model.at(earliestArray2, sequence[i]),
//                           model.at(earliestArray1, sequence[i]))),
//       model.iif(model.at(earliestArray2, sequence[i]) == -2,
//                 model.at(earliestArray1, sequence[i]),
//                 model.iif(prev + model.at(distanceArray, sequence[i - 1], sequence[i])
//                 >
//                               model.at(latestArray1, sequence[i]) -
//                                   model.at(serviceArray, sequence[i]),
//                           model.at(earliestArray2, sequence[i]),
//                           model.at(earliestArray1, sequence[i]))));
// });