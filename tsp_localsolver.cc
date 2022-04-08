

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
  LSExpression serviceSetUpDuration;
  LSExpression serviceQuantitiesMatrix;
  LSExpression serviceExclusionCost;
  LSExpression serviceMatrixIndex;

  LSExpression vehicleStartIndicies;
  LSExpression vehicleEndIndicies;
  LSExpression vehicleCapacitiesMatrix;

  map<string, int64> ids_map_;
  vector<map<int, LSExpression>> timesFromWarehouses;
  vector<map<int, LSExpression>> timesToWarehouses;

  LSExpression twStartsArray;
  LSExpression twEndsArray;
  LSExpression twAbsoluteEndsArray;
  LSExpression nbTwsArray;
  LSExpression waitNextTWArray;

  // LSExpression nextStart(const LSExpression& service, const LSExpression& time) {
  //   LSExpression timeSelector = model.createLambdaFunction([&](LSExpression tw_index) {
  //     return model.iif(model.at(twAbsoluteEndsArray, service, tw_index) >= time,
  //                      model.at(twStartsArray, service, tw_index),
  //                      model.at(twAbsoluteEndsArray, service, nbTwsArray[service] -
  //                      1));
  //   });

  //   LSExpression earliestAvailableTime =
  //       model.iif(nbTwsArray[service] == 0, 0,
  //                 model.min(model.range(0, nbTwsArray[service]), timeSelector));

  //   return model.max(time, earliestAvailableTime);
  // }

  LSExpression twEndSelect(const LSExpression& service, const LSExpression& time) {
    LSExpression timeWindowSelector =
        model.createLambdaFunction([&](LSExpression tw_index) {
          return model.iif(model.at(twAbsoluteEndsArray, service, tw_index) >= time &&
                               model.at(twStartsArray, service, tw_index) <= time,
                           model.at(twEndsArray, service, tw_index),
                           model.at(twEndsArray, service, nbTwsArray[service] - 1));
        });
    return model.min(model.range(0, nbTwsArray[service]), timeWindowSelector);
  }

  LSExpression nextStart(const LSExpression& service, const LSExpression& time) {
    LSExpression timeWindowSelector =
        model.createLambdaFunction([&](LSExpression tw_index) {
          LSExpression twDecisionAbsoluteEnd = model.iif(
              waitNextTWArray[service], model.at(twAbsoluteEndsArray, service, tw_index),
              model.at(twEndsArray, service, tw_index));
          return model.iif(
              twDecisionAbsoluteEnd >= time, model.at(twStartsArray, service, tw_index),
              model.at(twAbsoluteEndsArray, service, nbTwsArray[service] - 1));
        });
    LSExpression earliestAvailableTime =
        model.iif(nbTwsArray[service] == 0, 0,
                  model.min(model.range(0, nbTwsArray[service]), timeWindowSelector));
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
      , twAbsoluteEndsArray(model.array())
      , waitNextTWArray(model.array()) {
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
    vector<int> serviceSetUpDurationVec;
    vector<vector<float>> serviceQuantitiesVec;
    vector<float> serviceExclusionCostVec;
    vector<int> serviceMatrixIndexVec;
    vector<int> nbTWsVec;

    int s = 0;
    for (const auto& service : problem.services()) {
      serviceMatrixIndexVec.push_back(service.matrix_index());
      waitNextTWArray.addOperand(model.boolVar());
      serviceTimeVec.push_back(service.duration());
      serviceSetUpDurationVec.push_back(service.setup_duration());
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

      if (service.time_windows_size() > 0) {
        for (const auto& tw : service.time_windows()) {
          serviceTWStarts.push_back(tw.start());
          serviceTWEnds.push_back(tw.end());
          if (service.late_multiplier() > 0) {
            serviceTWAbsoluteEnds.push_back(tw.end() + tw.maximum_lateness());
          } else {
            serviceTWAbsoluteEnds.push_back(tw.end());
          }
        }
      } else {
        serviceTWStarts.push_back(0);
        serviceTWEnds.push_back(CUSTOM_MAX_INT);
        serviceTWAbsoluteEnds.push_back(CUSTOM_MAX_INT);
      }

      serviceTWAbsoluteEnds.back();
      twStartsArray.addOperand(
          model.array(serviceTWStarts.begin(), serviceTWStarts.end()));
      twStartsArray.setName("TW starts :");
      twEndsArray.addOperand(model.array(serviceTWEnds.begin(), serviceTWEnds.end()));
      twAbsoluteEndsArray.addOperand(
          model.array(serviceTWAbsoluteEnds.begin(), serviceTWAbsoluteEnds.end()));

      if (service.time_windows_size() == 0) {
        nbTWsVec.push_back(1);
      } else {
        nbTWsVec.push_back(service.time_windows_size());
      }
      ids_map_[(string)service.id()] = s;
      s++;
    }
    serviceMatrixIndex =
        model.array(serviceMatrixIndexVec.begin(), serviceMatrixIndexVec.end());
    serviceTime = model.array(serviceTimeVec.begin(), serviceTimeVec.end());
    serviceSetUpDuration =
        model.array(serviceSetUpDurationVec.begin(), serviceSetUpDurationVec.end());
    serviceExclusionCost =
        model.array(serviceExclusionCostVec.begin(), serviceExclusionCostVec.end());
    nbTwsArray = model.array(nbTWsVec.begin(), nbTWsVec.end());
  }

  //   void precedenceConstraints(LSModel& model) {}

  void createModelAndSolve() {
    // Decision Variables :

    // A tour for each vehicle
    vector<LSExpression> serviceSequences;
    serviceSequences.reserve(problem.vehicles_size() + 1);
    for (auto vehicle : problem.vehicles()) {
      serviceSequences.push_back(model.listVar(problem.services_size()));
      serviceSequences.back().setName("sequence_" + vehicle.id());
    }

    // extra_super_vehicle
    serviceSequences.push_back(model.listVar(problem.services_size()));
    serviceSequences.back().setName("sequence_unassigned");

    // Are the vehicles actually used?
    vector<LSExpression> vehiclesUsed;

    // Number of vehicles used
    LSExpression nbVehiclesUsed;
    vehiclesUsed.resize(problem.vehicles_size() + 1);

    vector<LSExpression> routeDuration(problem.vehicles_size());
    vector<LSExpression> endTime(problem.vehicles_size());
    vector<LSExpression> beginTime(problem.vehicles_size());
    vector<LSExpression> serviceStartsInTW(problem.vehicles_size());
    vector<LSExpression> lateness(problem.vehicles_size());
    vector<LSExpression> excessLateness(problem.vehicles_size());

    // all services must be satisfied by the vehicles
    LSExpression partition =
        model.partition(serviceSequences.begin(), serviceSequences.end());
    model.constraint(partition);

    LSExpression unassignedServices = serviceSequences[problem.vehicles_size()];
    LSExpression numberOfUnassignedServices = model.count(unassignedServices);
    vehiclesUsed[problem.vehicles_size()] = numberOfUnassignedServices > 0;

    // Constraints for each vehicle
    int k = 0;
    for (const auto& vehicle : problem.vehicles()) {
      LSExpression sequenceVehicle = serviceSequences[k];
      LSExpression c = model.count(sequenceVehicle);

      vehiclesUsed[k] = c > 0;
      vehiclesUsed[k].setName("vehicleUsed_" + to_string(k));

      LSExpression endSelector =
          model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
            return nextStart(
                       sequenceVehicle[i],
                       model.iif(
                           i == 0,
                           ((vehicle.has_time_window())
                                ? static_cast<lsint>(vehicle.time_window().start())
                                : 0) +
                               timesFromWarehouses[vehicle.matrix_index()]
                                                  [vehicle.start_index()]
                                                  [sequenceVehicle[i]] +
                               serviceSetUpDuration[sequenceVehicle[i]],
                           prev +
                               model.at(timeMatrices[vehicle.matrix_index()],
                                        sequenceVehicle[i - 1], sequenceVehicle[i]) +
                               model.iif(serviceMatrixIndex[sequenceVehicle[i]] ==
                                             serviceMatrixIndex[sequenceVehicle[i - 1]],
                                         0, serviceSetUpDuration[sequenceVehicle[i]]))) +
                   serviceTime[sequenceVehicle[i]];
          });

      endTime[k] = model.array(model.range(0, c), endSelector);

      LSExpression beginTimeSelector = model.createLambdaFunction([&](LSExpression i) {
        return endTime[k][i] - serviceTime[sequenceVehicle[i]];
        // to include setup duration of the first service in the begin time calculation
        // use the following one:
        // return endTime[k][i] - serviceTime[sequenceVehicle[i]] -
        //        model.iif(i > 0 && serviceMatrixIndex[sequenceVehicle[i]] ==
        //                               serviceMatrixIndex[sequenceVehicle[i - 1]],
        //                  0, serviceSetUpDuration[sequenceVehicle[i]]);
      });

      beginTime[k] = model.array(model.range(0, c), beginTimeSelector);

      routeDuration[k] =
          model.iif(c > 0,
                    endTime[k][c - 1] +
                        timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                         [sequenceVehicle[c - 1]],
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
        if (vehicle.cost_late_multiplier() > 0) {
          model.constraint(routeDuration[k] <=
                           static_cast<lsint>(vehicle.time_window().end() +
                                              vehicle.time_window().maximum_lateness()));
        } else {
          model.constraint(routeDuration[k] <=
                           static_cast<lsint>(vehicle.time_window().end()));
        }
      }
      LSExpression excessLatenessSelector =
          model.createLambdaFunction([&](LSExpression i) {
            // excess lateness is only calculated wrt the very last absolute end
            // because nextStart() can't return an infeasable intermediate time.
            return model.max(0, endTime[k][i] - serviceTime[sequenceVehicle[i]] -
                                    (model.at(twAbsoluteEndsArray, sequenceVehicle[i],
                                              nbTwsArray[sequenceVehicle[i]] - 1)));
          });
      excessLateness[k] = model.sum(model.range(0, c), excessLatenessSelector);

      LSExpression latenessSelector = model.createLambdaFunction([&](LSExpression i) {
        return model.max(
            0, endTime[k][i] - serviceTime[sequenceVehicle[i]] -
                   twEndSelect(sequenceVehicle[i],
                               endTime[k][i] - serviceTime[sequenceVehicle[i]]));
      });

      // latenessTW = model.range(0, nbTwsArray[sequenceVehicle[i]]);
      lateness[k] = model.sum(model.range(0, c), latenessSelector);

      for (int unit = 0; unit < vehicle.capacities_size(); unit++) {
        LSExpression quantityCumulator =
            model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
              return model.max(
                  0, prev + model.at(serviceQuantitiesMatrix, sequenceVehicle[i], unit));
            });
        LSExpression routeQuantityUnit =
            model.array(model.range(0, c), quantityCumulator);
        LSExpression quantityUnitChecker =
            model.createLambdaFunction([&](LSExpression i) {
              return routeQuantityUnit[i] <= vehicle.capacities(unit).limit();
            });
        model.constraint(model.and_(model.range(0, c), quantityUnitChecker));
      }

      for (const auto& relation : problem.relations()) {
        if (relation.type() == "shipment") {
          for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
               link_index++) {
            model.constraint(
                model.contains(sequenceVehicle, static_cast<lsint>(IdIndex(
                                                    relation.linked_ids(link_index)))) ==
                model.contains(
                    sequenceVehicle,
                    static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1)))));
            model.constraint(
                model.indexOf(sequenceVehicle, static_cast<lsint>(IdIndex(
                                                   relation.linked_ids(link_index)))) <=
                model.indexOf(
                    sequenceVehicle,
                    static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1)))));
          }
        }
        if (relation.type() == "order") {
          for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
               link_index++) {
            LSExpression sequenceContainsCurrentService = model.contains(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index))));
            LSExpression sequenceContainsNextService = model.contains(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1))));
            LSExpression nextIsNotAssigned = model.contains(
                unassignedServices,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1))));

            model.constraint(model.iif(sequenceContainsCurrentService,
                                       sequenceContainsNextService || nextIsNotAssigned,
                                       !sequenceContainsNextService));
            model.constraint(model.iif(
                sequenceContainsNextService,

                model.indexOf(sequenceVehicle, static_cast<lsint>(IdIndex(
                                                   relation.linked_ids(link_index)))) <=
                    model.indexOf(
                        sequenceVehicle,
                        static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1)))),
                true));
          }
        }
        if (relation.type() == "sequence") {
          for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
               link_index++) {
            LSExpression sequenceContainsCurrentService = model.contains(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index))));
            LSExpression sequenceContainsNextService = model.contains(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1))));
            LSExpression currentServiceIndexInSequence = model.indexOf(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index))));
            LSExpression nextServiceIndexInSequence = model.indexOf(
                sequenceVehicle,
                static_cast<lsint>(IdIndex(relation.linked_ids(link_index + 1))));

            model.constraint(sequenceContainsCurrentService ==
                             sequenceContainsNextService);
            model.constraint(model.iif(
                sequenceContainsCurrentService,
                currentServiceIndexInSequence + 1 == nextServiceIndexInSequence, true));
          }
        }
      }
      k++;
    }

    LSExpression totalExcessLateness =
        model.sum(excessLateness.begin(), excessLateness.end());
    model.constraint(totalExcessLateness == 0);

    LSExpression totalLateness = model.sum(lateness.begin(), lateness.end());
    totalLateness.setName("total Lateness");

    LSExpression exclusionCostCumulator = model.createLambdaFunction(
        [&](LSExpression i) { return serviceExclusionCost[unassignedServices[i]]; });

    LSExpression totalExclusionCost =
        model.sum(model.range(0, numberOfUnassignedServices), exclusionCostCumulator);
    totalExclusionCost.setName("total Exclusion Cost");

    // Total vehicles used :
    nbVehiclesUsed = model.sum(vehiclesUsed.begin(), vehiclesUsed.end());

    // Total distance traveled :
    LSExpression totalDuration;
    totalDuration = model.sum(routeDuration.begin(), routeDuration.end());
    totalDuration.setName("totalDuration");

    // model.minimize(nbVehiclesUsed);
    // model.minimize(totalExcessLateness);
    model.minimize(totalExclusionCost);
    model.minimize(totalLateness);
    model.minimize(totalDuration);

    model.close();

    cout << model.toString() << endl;
    localsolver.getParam().setTimeLimit(FLAGS_time_limit_in_ms / 1000);

    // With ls a LocalSolver object
    auto iis = localsolver.computeInconsistency();
    cout << iis.toString() << endl;

    localsolver.solve();

    cout << "-------------------- DATA ------------------------------" << endl;

    for (auto vehicle : problem.vehicles()) {
      cout << "timeMatrix of : " << vehicle.id() << " "
           << timeMatrices[vehicle.matrix_index()].getArrayValue().toString() << endl;
      cout << "timeFromWarehouses of : " << vehicle.id() << " "
           << timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                  .getArrayValue()
                  .toString()
           << endl;
      cout << "timeToWarehouses of : " << vehicle.id() << " "
           << timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                  .getArrayValue()
                  .toString()
           << endl;

      cout << "setup duration of services :"
           << serviceSetUpDuration.getArrayValue().toString() << endl;
      cout << "maximum duration of " << vehicle.id() << " : " << vehicle.duration()
           << endl;
    }

    cout << "number of time windows " << nbTwsArray.getArrayValue().toString() << endl;

    cout << " Time Window Starts Array :  " << twStartsArray.getArrayValue().toString()
         << endl;

    cout << " Time Window Ends Array :  " << twEndsArray.getArrayValue().toString()
         << endl;
    cout << " Time Window Absolute End Array :  "
         << twAbsoluteEndsArray.getArrayValue().toString() << endl;

    cout << endl;
    cout << "total duration = " << totalDuration.getValue() << endl;

    cout << "-----------------------------------------------------------------" << endl;

    cout << "vehicle capacities matrix"
         << vehicleCapacitiesMatrix.getArrayValue().toString() << endl;

    cout << "service Quantities Matrix "
         << serviceQuantitiesMatrix.getArrayValue().toString() << endl;
    cout << " ====================== SOLUTION ==================================="
         << endl;
    cout << "--------------------   ASSIGNEMENTS  ----------------------------" << endl;

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
    }
    cout << " ----------------------begin times--------------------------  " << endl;
    for (int k = 0; k < problem.vehicles_size(); k++) {
      cout << "begin times service of " << problem.vehicles(k).id() << " "
           << beginTime[k].getArrayValue().toString() << endl;
    }
    cout << " ----------------------End times--------------------------  " << endl;
    for (int k = 0; k < problem.vehicles_size(); k++) {
      cout << "end times service of " << problem.vehicles(k).id() << " "
           << endTime[k].getArrayValue().toString() << endl;
      cout << "total route duration :  " << problem.vehicles(k).id() << " "
           << routeDuration[k].getValue() << endl;
    }
    cout << endl;
    cout << " -------------------- LATENESS -----------------------------------" << endl;
    cout << " total Absolute Latenness " << totalExcessLateness.getValue() << endl;
    cout << " total Latenness " << totalLateness.getValue() << endl;
    cout << endl;

    cout << " -------------------- OBJECTIVE -----------------------------------" << endl;
    for (int i = 0; i < model.getNbObjectives(); i++) {
      cout << model.getObjective(i).toString() << endl;
    }
    cout << endl;
    cout << " ----------------------- END -------------------------------------------"
         << endl;
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