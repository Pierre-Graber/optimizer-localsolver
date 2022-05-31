

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <lsarray.h>
#include <lscallbacktype.h>
#include <lscollection.h>
#include <lsexpression.h>
#include <lssolution.h>
#include <ostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

// https://stackoverflow.com/a/77336
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "./localsolver_result.pb.h"
#include "./localsolver_vrp.pb.h"
// #include "callback_ls.cpp"
#include "localsolver.h"
#include "ortools/base/commandlineflags.h"

#define CUSTOM_MAX_INT (int64) std::pow(2, 30)

#define CUSTOM_BIGNUM_COST 1e6

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

void handler(int sig) {
  void* array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

bool equals(double a, double b) {
  float EPSILON = 10e-6;
  return fabs(a - b) < EPSILON;
}

class MyCallback : public LSCallback {
private:
  int lastBestRunningTime;
  lsdouble lastBestValue;
  const int64 minimum_duration_in_s;
  const int64 time_out_multiplier;

public:
  MyCallback(const int64 minimum_duration, const int64 time_out_multip)
      : minimum_duration_in_s(minimum_duration / 1000)
      , time_out_multiplier(time_out_multip) {
    lastBestRunningTime = 0;
    lastBestValue = CUSTOM_MAX_INT;
  }

  void callback(LocalSolver& ls, LSCallbackType type) {
    (void)type;
    LSStatistics stats = ls.getStatistics();
    LSExpression obj = ls.getModel().getObjective(0);

    if (obj.getDoubleValue() < lastBestValue) {
      lastBestRunningTime = stats.getRunningTime();
      lastBestValue = obj.getDoubleValue();
    }
    cout << "getRunningTime " << stats.getRunningTime() << " " << lastBestRunningTime
         << " " << lastBestValue << endl;

    if (stats.getRunningTime() > minimum_duration_in_s &&
        stats.getRunningTime() > time_out_multiplier * lastBestRunningTime) {
      cout << ">>>>>>>  Resolution stopped" << endl;
      ls.stop();
    }
  }
};
class localsolver_VRP {
public:
  const localsolver_vrp::Problem problem;

  LocalSolver localsolver;
  LSModel model;
  vector<LSExpression> timeMatrices; // the matrices are reordered by service.matrix_index
  vector<LSExpression>
      distanceMatrices; // the matrices are reordered by service.matrix_index

  LSExpression serviceTime;
  LSExpression serviceLateMultiplier;
  LSExpression serviceSetUpDuration;
  LSExpression serviceQuantitiesMatrix;
  LSExpression serviceExclusionCost;
  LSExpression servicePriority;
  LSExpression serviceMatrixIndex;
  int maxTwStarts;

  LSExpression vehicleStartIndicies;
  LSExpression vehicleEndIndicies;
  LSExpression vehicleCapacitiesMatrix;

  LSExpression endOfVehicleTimeWindow;

  map<string, int64> service_ids_map_;
  vector<map<int, LSExpression>> timesFromWarehouses;
  vector<map<int, LSExpression>> timesToWarehouses;

  vector<map<int, LSExpression>> distanceFromWarehouses;
  vector<map<int, LSExpression>> distanceToWarehouses;

  LSExpression serviceTwStartsArray;
  LSExpression serviceTwEndsArray;
  LSExpression twAbsoluteEndsArray;
  LSExpression nbTwsArray;
  LSExpression waitNextTWArray;

  vector<LSExpression> startTimeVehicle;
  vector<LSExpression> routeDuration;
  vector<LSExpression> routeDurationCost;
  vector<LSExpression> routeDistance;
  vector<LSExpression> routeDistanceCost;
  vector<LSExpression> routeFixedCost;
  vector<LSExpression> endTime;
  vector<LSExpression> endTimeWithRest;
  vector<LSExpression> beginTime;
  vector<LSExpression> arrivalTime;
  vector<LSExpression> waitingTime;
  vector<LSExpression> serviceStartsInTW;
  vector<LSExpression> latenessCost;
  vector<LSExpression> latenessOfServicesOfVehicle;
  vector<LSExpression> excessLateness;
  vector<LSExpression> Rest;
  vector<LSExpression> restDuration;
  vector<int> allVehicleIndices;
  map<string, int64> vehicle_ids_map_;
  vector<map<string, int64>> rest_ids_map;

  // Time Leaving The Warehouse
  vector<LSExpression> timeLeavingTheWarehouse;

  LSExpression twEndSelect(const LSExpression& service, const LSExpression& time) {
    LSExpression timeWindowSelector =
        model.createLambdaFunction([&](LSExpression tw_index) {
          return model.iif(
              model.at(twAbsoluteEndsArray, service, tw_index) >= time &&
                  model.at(serviceTwStartsArray, service, tw_index) <= time &&
                               waitNextTWArray[service] == 0,
              model.at(serviceTwEndsArray, service, tw_index),
              model.at(serviceTwEndsArray, service, nbTwsArray[service] - 1));
        });
    return model.min(model.range(0, nbTwsArray[service]), timeWindowSelector);
  }

  LSExpression nextStart(const LSExpression& service, const LSExpression& time) {
    LSExpression timeWindowSelector =
        model.createLambdaFunction([&](LSExpression tw_index) {
          LSExpression twDecisionAbsoluteEnd =
              model.iif(waitNextTWArray[service] == 1,
                        model.at(serviceTwEndsArray, service, tw_index),
              model.at(twAbsoluteEndsArray, service, tw_index));
          return model.iif(
              twDecisionAbsoluteEnd >= time,
              model.at(serviceTwStartsArray, service, tw_index),
              model.at(twAbsoluteEndsArray, service, nbTwsArray[service] - 1));
        });
    LSExpression earliestAvailableTime =
        model.iif(nbTwsArray[service] == 0, 0,
                  model.min(model.range(0, nbTwsArray[service]), timeWindowSelector));
    return model.max(time, earliestAvailableTime);
  }

  void firstAndSecondSolving(vector<LSExpression> timeLeavingTheWarehouseConstraint) {
    cout << model.toString() << endl;
    std::unique_ptr<MyCallback> cbp;
    if (!FLAGS_only_first_solution) {
      cbp =
          std::make_unique<MyCallback>(FLAGS_minimum_duration, FLAGS_time_out_multiplier);
      localsolver.addCallback(localsolver::CT_TimeTicked, cbp.get());
      localsolver.getParam().setTimeLimit((FLAGS_time_limit_in_ms / 1000) / 2);

      auto iis = localsolver.computeInconsistency();
      cout << iis.toString() << endl;

      localsolver.solve();
    } else {
      // localsolver.getParam().setIterationLimit(1);

      auto iis = localsolver.computeInconsistency();
      cout << iis.toString() << endl;
      cbp = std::make_unique<MyCallback>(3000, FLAGS_time_out_multiplier);
      localsolver.addCallback(localsolver::CT_TimeTicked, cbp.get());
    }

    model.open();
    for (auto constraint : timeLeavingTheWarehouseConstraint) {
      model.removeConstraint(constraint);
    }
    model.close();

    localsolver.getParam().setTimeLimit((FLAGS_time_limit_in_ms / 1000) / 2);

    localsolver.solve();
  }

  void RelationBuilder(LSExpression& sequenceVehicle, LSExpression& unassignedServices) {
    for (const auto& relation : problem.relations()) {
      if (relation.type() == "shipment") {
        for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
             link_index++) {
          model.constraint(
              model.contains(sequenceVehicle,
                             static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                        service_ids_map_))) ==
              model.contains(sequenceVehicle, static_cast<lsint>(IdIndex(
                                                  relation.linked_ids(link_index + 1),
                                                  service_ids_map_))));
          model.constraint(
              model.indexOf(sequenceVehicle,
                            static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                       service_ids_map_))) <=
              model.indexOf(sequenceVehicle,
                            static_cast<lsint>(IdIndex(
                                relation.linked_ids(link_index + 1), service_ids_map_))));
        }
      }
      if (relation.type() == "order") {
        for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
             link_index++) {
          LSExpression sequenceContainsCurrentService = model.contains(
              sequenceVehicle, static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                          service_ids_map_)));
          LSExpression sequenceContainsNextService =
              model.contains(sequenceVehicle,
                             static_cast<lsint>(IdIndex(
                                 relation.linked_ids(link_index + 1), service_ids_map_)));
          LSExpression nextIsNotAssigned =
              model.contains(unassignedServices,
                             static_cast<lsint>(IdIndex(
                                 relation.linked_ids(link_index + 1), service_ids_map_)));

          model.constraint(model.iif(sequenceContainsCurrentService,
                                     sequenceContainsNextService || nextIsNotAssigned,
                                     !sequenceContainsNextService));
          model.constraint(model.iif(
              sequenceContainsNextService,

              model.indexOf(sequenceVehicle,
                            static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                       service_ids_map_))) <=
                  model.indexOf(sequenceVehicle, static_cast<lsint>(IdIndex(
                                                     relation.linked_ids(link_index + 1),
                                                     service_ids_map_))),
              true));
        }
      }
      if (relation.type() == "sequence") {
        for (int link_index = 0; link_index < relation.linked_ids_size() - 1;
             link_index++) {
          LSExpression sequenceContainsCurrentService = model.contains(
              sequenceVehicle, static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                          service_ids_map_)));
          LSExpression sequenceContainsNextService =
              model.contains(sequenceVehicle,
                             static_cast<lsint>(IdIndex(
                                 relation.linked_ids(link_index + 1), service_ids_map_)));
          LSExpression indexOfCurrentServiceInSequence = model.indexOf(
              sequenceVehicle, static_cast<lsint>(IdIndex(relation.linked_ids(link_index),
                                                          service_ids_map_)));
          LSExpression indexOfNextServiceInSequence =
              model.indexOf(sequenceVehicle,
                            static_cast<lsint>(IdIndex(
                                relation.linked_ids(link_index + 1), service_ids_map_)));

          model.constraint(model.iif(
              sequenceContainsNextService,
              sequenceContainsCurrentService == sequenceContainsNextService, true));
          model.constraint(model.iif(
              sequenceContainsNextService,
              indexOfCurrentServiceInSequence + 1 == indexOfNextServiceInSequence, true));
        }
      }
    }
  }
  void capacityConstraintsOfVehicle(int k, LSExpression sequenceVehicle, LSExpression c) {
    const localsolver_vrp::Vehicle& vehicle = problem.vehicles(k);
    for (int unit = 0; unit < vehicle.capacities_size(); unit++) {
      LSExpression quantityCumulator =
          model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
            return model.max(
                0, prev + model.at(serviceQuantitiesMatrix, sequenceVehicle[i], unit));
          });
      LSExpression routeQuantityUnit = model.array(model.range(0, c), quantityCumulator);
      LSExpression quantityUnitChecker = model.createLambdaFunction([&](LSExpression i) {
        return routeQuantityUnit[i] <= vehicle.capacities(unit).limit();
      });
      model.constraint(model.and_(model.range(0, c), quantityUnitChecker));
    }
  }

  void MatrixBuilder(vector<LSExpression>& Matrices, const RepeatedField<float>& matrix) {
    LSExpression Matrix(model.array());
    int matrix_size = sqrt(matrix.size());

    for (const auto& serviceFrom : problem.services()) {
      vector<int> matrix_row;
      matrix_row.reserve(matrix_size);
      for (const auto& serviceTo : problem.services()) {
        if (matrix.size() > 0) {
          matrix_row.push_back(ceil(matrix.at(serviceFrom.matrix_index() * matrix_size +
                                              serviceTo.matrix_index())));
        } else {
          matrix_row.push_back(0);
        }
      }
      Matrix.addOperands(model.array(matrix_row.begin(), matrix_row.end()));
    }
    Matrices.push_back(Matrix);
  }

  int64 IdIndex(const string& id, map<string, int64> map_) const {
    map<string, int64>::const_iterator it = map_.find(id);
    if (it != map_.end())
      return it->second;
    else
      return -1;
  }

  string IndexId(const int index, map<string, int64> map_) const {
    for (auto it = map_.begin(); it != map_.end(); ++it) {
      if (it->second == index) {
        return it->first;
      }
    }
    return "not found";
  }

  LSExpression needsPause(LSExpression restStart, LSExpression start, LSExpression end) {
    return start <= restStart && end >= restStart;
  }

  LSExpression travelAndSetupEndWithPause(localsolver_vrp::Vehicle vehicle,
                                          LSExpression i, LSExpression time,
                                          vector<LSExpression> serviceSequences) {
    int k = IdIndex(vehicle.id(), vehicle_ids_map_);
    LSExpression sequenceVehicle = serviceSequences[k];
    LSExpression travelAndSetupDuration =
        model.iif(i == 0,
                  timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                     [sequenceVehicle[0]] +
                      serviceSetUpDuration[sequenceVehicle[0]],
                  model.at(timeMatrices[vehicle.matrix_index()], sequenceVehicle[i - 1],
                           sequenceVehicle[i]) +
                      model.iif(serviceMatrixIndex[sequenceVehicle[i]] ==
                                    serviceMatrixIndex[sequenceVehicle[i - 1]],
                                0, serviceSetUpDuration[i]));
    LSExpression travelAndSetupEnd = time + travelAndSetupDuration;
    LSExpression extraPauseDurationSelector =
        model.createLambdaFunction([&](LSExpression rest_index) {
          return model.iif(needsPause(Rest[k][rest_index], time, travelAndSetupEnd),
                           restDuration[k][rest_index], 0);
        });
    LSExpression extraPauseDuration =
        model.sum(model.range(0, vehicle.rests_size()), extraPauseDurationSelector);
    return travelAndSetupEnd + extraPauseDuration;
  }

  LSExpression serviceEndWithPause(localsolver_vrp::Vehicle vehicle, LSExpression service,
                                   LSExpression time) {
    LSExpression nextStartWithoutPause = nextStart(service, time);
    LSExpression endWithoutPause = nextStartWithoutPause + serviceTime[service];
    LSExpression firstPauseselector =
        model.createLambdaFunction([&](LSExpression rest_index) {
          return model.iif(
              needsPause(
                  model.at(Rest[IdIndex(vehicle.id(), vehicle_ids_map_)], rest_index),
                  time, endWithoutPause),
              model.at(Rest[IdIndex(vehicle.id(), vehicle_ids_map_)], rest_index),
              static_cast<lsint>(CUSTOM_MAX_INT));
        });
    LSExpression restIndexSelector =
        model.createLambdaFunction([&](LSExpression rest_index) {
          return model.iif(
              needsPause(
                  model.at(Rest[IdIndex(vehicle.id(), vehicle_ids_map_)], rest_index),
                  time, endWithoutPause),
              rest_index, 0);
        });
    LSExpression firstPause =
        model.min(model.range(0, vehicle.rests_size()), firstPauseselector);
    LSExpression rest_index =
        model.max(model.range(0, vehicle.rests_size()), restIndexSelector);
    LSExpression nextStartWithPause = nextStart(
        service,
        firstPause + restDuration[IdIndex(vehicle.id(), vehicle_ids_map_)][rest_index]);
    LSExpression endWithPause = nextStartWithPause + serviceTime[service];
    return model.iif(firstPause == static_cast<lsint>(CUSTOM_MAX_INT), endWithoutPause,
                     endWithPause);
  }

  void setInitialSolution(vector<LSExpression>& serviceSequences) {
    if (FLAGS_only_first_solution) {
      for (auto const& route : problem.routes()) {
        int vehicle_index = IdIndex(route.vehicle_id(), vehicle_ids_map_);
        LSExpression sequenceVehicle = serviceSequences[vehicle_index];
        // cout << "service_id :" << route.service_ids(0) << endl;
        for (int index_of_service_in_route = 0;
             index_of_service_in_route < route.service_ids_size() - 1;
             index_of_service_in_route++) {
          // model.constraint(model.contains(sequenceVehicle, service));
          // model.constraint(model.indexOf(sequenceVehicle, service) ==
          //                  index_of_service_in_route);

          LSExpression sequenceContainsCurrentService = model.contains(
              sequenceVehicle,
              static_cast<lsint>(IdIndex(route.service_ids(index_of_service_in_route),
                                         service_ids_map_)));
          LSExpression sequenceContainsNextService = model.contains(
              sequenceVehicle,
              static_cast<lsint>(IdIndex(route.service_ids(index_of_service_in_route + 1),
                                         service_ids_map_)));
          LSExpression indexOfCurrentServiceInSequence = model.indexOf(
              sequenceVehicle,
              static_cast<lsint>(IdIndex(route.service_ids(index_of_service_in_route),
                                         service_ids_map_)));
          LSExpression indexOfNextServiceInSequence = model.indexOf(
              sequenceVehicle,
              static_cast<lsint>(IdIndex(route.service_ids(index_of_service_in_route + 1),
                                         service_ids_map_)));
          model.constraint(sequenceContainsCurrentService == sequenceContainsNextService);
          model.constraint(indexOfCurrentServiceInSequence == index_of_service_in_route);
          // model.constraint(model.iif(
          //     sequenceContainsCurrentService,
          //     indexOfCurrentServiceInSequence + 1 == indexOfNextServiceInSequence,
          //     true));

          // model.constraint(sequenceContainsCurrentService);
          // model.constraint(sequenceContainsNextService);
          // model.constraint(indexOfCurrentServiceInSequence + 1 ==
          //                  indexOfNextServiceInSequence);
        }
      }
    }

    model.close();
    LSSolution sol = localsolver.getSolution();
    std::unordered_set<string> initializedServiceIds;
    map<string, localsolver_vrp::TimeWindow> used_tw_for_service_map;
    for (const auto& route : problem.routes()) {
      int route_index = IdIndex(route.vehicle_id(), vehicle_ids_map_);
      const localsolver_vrp::Vehicle& vehicle = problem.vehicles(route_index);
      if (vehicle.shift_preference() == "force_start") {
        LSExpression tLTW = localsolver.getModel().getExpression(
            "timeLeavingTheWarehouse" + vehicle.id());
        sol.setValue(tLTW, static_cast<lsint>(vehicle.time_window().start()));
      } else {
        vector<int> start_time;
        start_time.reserve(route.service_ids().size());
        const RepeatedField timeMatrix = problem.matrices(vehicle.matrix_index()).time();
        LSExpression listExpr =
            localsolver.getModel().getExpression("sequence_" + vehicle.id());
        LSCollection sequence = listExpr.getCollectionValue();

        int time_matrix_size = sqrt(problem.matrices(vehicle.matrix_index()).time_size());
        int previous_end = static_cast<int>(vehicle.time_window().start()) | 0;
        int previous_location_index = vehicle.start_index();
        for (const auto& service_id : route.service_ids()) {
          const localsolver_vrp::Service& service =
              problem.services(IdIndex(service_id, service_ids_map_));
          uint current_arrival =
              previous_end +
              timeMatrix.at(previous_location_index * time_matrix_size +
                            service.matrix_index()) +
              service.setup_duration();
          int current_start = current_arrival;
          int current_location_index = service.matrix_index();
          int tw_index = 0;
          if (service.time_windows_size() == 0) {
            current_start = current_arrival;
          } else {
            for (const auto& tw : service.time_windows()) {
              if (current_arrival <= (service.late_multiplier() > 0
                                          ? tw.end() + tw.maximum_lateness()
                                          : tw.end())) {
                current_start = max<uint>(current_arrival, tw.start());
                used_tw_for_service_map[service_id] = tw;
                break;
              }
              tw_index++;
            }
          }
          start_time.push_back(current_start);
          previous_end = current_start + service.duration();
          previous_location_index = current_location_index;
          if (service_ids_map_.count(service_id)) {
            sequence.add(IdIndex(service_id, service_ids_map_));
            initializedServiceIds.insert(service_id);
          }
        }

        localsolver_vrp::Service next_service = problem.services(
            IdIndex(route.service_ids(start_time.size() - 1), service_ids_map_));
        localsolver_vrp::Service current_service;
        for (int service_index = start_time.size() - 2; service_index >= 0;
             service_index--) {
          current_service = problem.services(
              IdIndex(route.service_ids(service_index), service_ids_map_));
          int time_between_two_starts =
              start_time[service_index + 1] - start_time[service_index];

          int idle_time =
              time_between_two_starts -
              (timeMatrix.at(current_service.matrix_index() * time_matrix_size +
                             next_service.matrix_index()) +
               current_service.duration() +
               (current_service.matrix_index() == next_service.matrix_index()
                    ? 0
                    : next_service.setup_duration()));
          const localsolver_vrp::TimeWindow& tw_used =
              used_tw_for_service_map.find(current_service.id())->second;
          if (idle_time > 0) {
            if (used_tw_for_service_map.find(current_service.id()) ==
                used_tw_for_service_map.end()) {
              start_time[service_index] += idle_time;
            } else {
              // cout << current_service.id()
              //      << " initial start time : " << start_time[service_index] << endl;
              start_time[service_index] +=
                  min<int>(idle_time, tw_used.end() +
                                          (current_service.late_multiplier() > 0
                                               ? tw_used.maximum_lateness()
                                               : 0) -
                                          start_time[service_index]);
              cout << " mooved start time : " << start_time[service_index] << endl;
            }
          }
          next_service = current_service;
        }
        int setTimeLeavingWarehouse =
            start_time[0] -
            timeMatrix.at(
                vehicle.start_index() * time_matrix_size +
                problem.services(IdIndex(route.service_ids(0), service_ids_map_))
                    .matrix_index()) -
            current_service.setup_duration();
        cout << "setTimeLeavingWarehouse : " << setTimeLeavingWarehouse << endl;
        ;
        LSExpression tLTW = localsolver.getModel().getExpression(
            "timeLeavingTheWarehouse" + vehicle.id());
        sol.setValue(tLTW, static_cast<lsint>(setTimeLeavingWarehouse));
      }
    }
    LSCollection sequenceUnassigned =
        localsolver.getModel().getExpression("sequence_unassigned").getCollectionValue();

    for (const auto& service : problem.services()) {
      if (initializedServiceIds.count(service.id()) == 0)
        sequenceUnassigned.add(IdIndex(service.id(), service_ids_map_));
    }
  }

  void ParseSolution(localsolver_result::Result* result,
                     const vector<LSExpression>& servicesSequence,
                     const vector<LSExpression> vehicleUsed) {
    LSStatistics stats = localsolver.getStatistics();
    long nbOfIterations = stats.getNbIterations();
    int runningTime = stats.getRunningTime();
    result->clear_routes();

    for (int route_index = 0; route_index < problem.vehicles_size(); route_index++) {
      if (vehicleUsed[route_index].getValue() == 1) {
      LSArray beginTimeArray = beginTime[route_index].getArrayValue();
      LSArray endTimeArray = endTime[route_index].getArrayValue();
      LSArray latenessServicesOfVehicleArray =
          latenessOfServicesOfVehicle[route_index].getArrayValue();
      LSArray timeToWareHouseArray =
          timesToWarehouses[problem.vehicles(route_index).matrix_index()]
                           [problem.vehicles(route_index).end_index()]
                               .getArrayValue();
      LSCollection servicesCollection =
          servicesSequence[route_index].getCollectionValue();
      LSArray restBeginTimeArray = Rest[route_index].getArrayValue();
      LSArray restDurationArray = restDuration[route_index].getArrayValue();
        localsolver_result::Route* route = result->add_routes();
        if (problem.vehicles(route_index).start_index() != -1) {
          localsolver_result::Activity* start_route = route->add_activities();
          start_route->set_type("start");
          start_route->set_index(-1);
          start_route->set_start_time(startTimeVehicle[route_index].getIntValue());
        }
        int currentRestIndex = 0;
        for (int activity_index = 0; activity_index < servicesCollection.count();
             activity_index++) {
          if (problem.vehicles(route_index).rests_size() > 0) {
            while (restBeginTimeArray.count() > currentRestIndex &&
                   beginTimeArray.getIntValue(activity_index) >=
                       restBeginTimeArray.getIntValue(currentRestIndex)) {
              localsolver_result::Activity* rest = route->add_activities();
              rest->set_type("break");
              rest->set_start_time(restBeginTimeArray.getIntValue(currentRestIndex));
              rest->set_id(IndexId(currentRestIndex, rest_ids_map[route_index]));
              currentRestIndex++;
            }
          }
          localsolver_result::Activity* activity = route->add_activities();
          activity->set_type("service");
          activity->set_id(IndexId(servicesCollection[activity_index], service_ids_map_));
          activity->set_index(
              problem.services(servicesCollection[activity_index]).problem_index());
          activity->set_start_time(beginTimeArray.getIntValue(activity_index));
          activity->set_lateness(
              latenessServicesOfVehicleArray.getDoubleValue(activity_index));
          for (auto const& quantity :
               problem.services(servicesCollection[activity_index]).quantities()) {
            activity->add_quantities(quantity);
          }
        }
        if (problem.vehicles(route_index).rests_size() > 0) {
          if (restBeginTimeArray.count() > currentRestIndex) {
            while (currentRestIndex < restBeginTimeArray.count()) {
              localsolver_result::Activity* rest = route->add_activities();
              rest->set_type("break");
              rest->set_start_time(restBeginTimeArray.getIntValue(currentRestIndex));
              rest->set_id(IndexId(currentRestIndex, rest_ids_map[route_index]));
              currentRestIndex++;
            }
          }
        }
        auto route_costs = route->mutable_cost_details();
        route_costs->set_fixed(problem.vehicles(route_index).cost_fixed());
        route_costs->set_distance(routeDistanceCost[route_index].getDoubleValue());
        route_costs->set_time(routeDurationCost[route_index].getDoubleValue());
        route_costs->set_lateness(latenessCost[route_index].getDoubleValue());
        if (problem.vehicles(route_index).end_index() != -1) {
          localsolver_result::Activity* end_route = route->add_activities();
          end_route->set_type("end");
          end_route->set_index(-1);
          end_route->set_start_time(
              max(endTimeArray.getIntValue(servicesCollection.count() - 1),
                  restBeginTimeArray.getIntValue(restBeginTimeArray.count() - 1) +
                      restDurationArray.getIntValue(restDurationArray.count() - 1)) +
              timeToWareHouseArray.getIntValue(
                  servicesCollection[servicesCollection.count() - 1]));
        }
      } else {
        localsolver_result::Route* route = result->add_routes();
        // if (problem.vehicles(route_index).start_index() != -1) {
          localsolver_result::Activity* start_route = route->add_activities();
          start_route->set_type("start");
          start_route->set_index(-1);
          start_route->set_start_time(startTimeVehicle[route_index].getIntValue());
        // }
        // if (problem.vehicles(route_index).end_index() != -1) {
          localsolver_result::Activity* end_route = route->add_activities();
          end_route->set_type("end");
          end_route->set_index(-1);
        end_route->set_start_time(endTimeVehicle[route_index].getIntValue());
        // }
      }
    }
    result->set_duration(runningTime);
    result->set_iterations(nbOfIterations);
  }

  // Constructor
  localsolver_VRP(localsolver_vrp::Problem& pb)
      : problem(pb)
      , localsolver()
      , model(localsolver.getModel())
      , timeMatrices(0, model.array())
      , distanceMatrices(0, model.array())
      , serviceQuantitiesMatrix(model.array())
      , maxTwStarts(0)
      , vehicleCapacitiesMatrix(model.array())
      , endOfVehicleTimeWindow(model.array())
      , serviceTwStartsArray(model.array())
      , serviceTwEndsArray(model.array())
      , twAbsoluteEndsArray(model.array())
      , waitNextTWArray(model.array())
      , startTimeVehicle(problem.vehicles_size())
      , routeDuration(problem.vehicles_size())
      , routeDurationCost(problem.vehicles_size())
      , routeDistance(problem.vehicles_size())
      , routeDistanceCost(problem.vehicles_size())
      , routeFixedCost(problem.vehicles_size())
      , endTime(problem.vehicles_size())
      , endTimeWithRest(problem.vehicles_size())
      , beginTime(problem.vehicles_size())
      , arrivalTime(problem.vehicles_size())
      , waitingTime(problem.vehicles_size())
      , serviceStartsInTW(problem.vehicles_size())
      , latenessCost(problem.vehicles_size())
      , latenessOfServicesOfVehicle(problem.vehicles_size())
      , excessLateness(problem.vehicles_size())
      , Rest(problem.vehicles_size(), model.array())
      , restDuration(problem.vehicles_size(), model.array())
      , rest_ids_map(problem.vehicles_size())
      , timeLeavingTheWarehouse(problem.vehicles_size()) {
    for (const auto& matrix : problem.matrices()) {
      MatrixBuilder(timeMatrices, matrix.time());
      MatrixBuilder(distanceMatrices, matrix.distance());
    }

    vector<int> vehicleStartIndiciesVec;
    vector<int> vehicleEndIndiciesVec;
    vector<vector<float>> vehicleCapacitiesMatrixVec;
    int k = 0;
    timesFromWarehouses.resize(problem.matrices_size());
    timesToWarehouses.resize(problem.matrices_size());
    distanceFromWarehouses.resize(problem.matrices_size());
    distanceToWarehouses.resize(problem.matrices_size());
    for (const auto& vehicle : problem.vehicles()) {
      endOfVehicleTimeWindow.addOperand(
          vehicle.has_time_window() ? static_cast<lsint>(vehicle.time_window().end())
                                    : CUSTOM_MAX_INT);
      vector<int> restDurationVec;
      vector<LSExpression> restVec;
      if (!vehicle.rests().empty()) {
        int rest_index = 0;
        for (auto const& rest : vehicle.rests()) {
          restVec.push_back(
              model.intVar(rest.time_window().start(), rest.time_window().end()));
          restDurationVec.push_back(rest.duration());
          rest_ids_map[k][(string)rest.id()] = rest_index;
          rest_index++;
        }
      } else {
        restVec.push_back(
            model.intVar(vehicle.time_window().start(), vehicle.time_window().start()));
        restDurationVec.push_back(0);
      }
      restDuration[k] = model.array(restDurationVec.begin(), restDurationVec.end());
      Rest[k] = model.array(restVec.begin(), restVec.end());

      int time_matrix_size = sqrt(problem.matrices(vehicle.matrix_index()).time_size());
      int distance_matrix_size =
          sqrt(problem.matrices(vehicle.matrix_index()).distance_size());
      map<int, LSExpression>& vehicleStartMapTime =
          timesFromWarehouses[vehicle.matrix_index()];
      map<int, LSExpression>& vehicleStartMapDistance =
          distanceFromWarehouses[vehicle.matrix_index()];
      if (vehicleStartMapTime.find(vehicle.start_index()) == vehicleStartMapTime.end()) {
        vector<int> timesFromWarehouse;
        timesFromWarehouse.reserve(problem.services_size());
        if (vehicle.start_index() > -1) {
          if (problem.matrices(vehicle.matrix_index()).time().size() > 0) {
            for (const auto& service : problem.services()) {
              timesFromWarehouse.push_back(
                  problem.matrices(vehicle.matrix_index())
                      .time(vehicle.start_index() * time_matrix_size +
                            service.matrix_index()));
            }
          } else {
            for (const auto& service : problem.services()) {
              (void)service;
              timesFromWarehouse.push_back(0);
            }
          }
        } else {
          timesFromWarehouse.resize(problem.services_size());
          fill_n(timesFromWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleStartMapTime[vehicle.start_index()] =
            model.array(timesFromWarehouse.begin(), timesFromWarehouse.end());
      }
      if (vehicleStartMapDistance.find(vehicle.start_index()) ==
          vehicleStartMapDistance.end()) {
        vector<int> distanceFromWarehouse;
        distanceFromWarehouse.reserve(problem.services_size());
        if (vehicle.start_index() > -1) {
          if (problem.matrices(vehicle.matrix_index()).distance().size() > 0) {
            for (const auto& service : problem.services()) {
              distanceFromWarehouse.push_back(
                  problem.matrices(vehicle.matrix_index())
                      .distance(vehicle.start_index() * distance_matrix_size +
                                service.matrix_index()));
            }
          } else {
            for (const auto& service : problem.services()) {
              (void)service;
              distanceFromWarehouse.push_back(0);
            }
          }
        } else {
          distanceFromWarehouse.resize(problem.services_size());
          fill_n(distanceFromWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleStartMapDistance[vehicle.start_index()] =
            model.array(distanceFromWarehouse.begin(), distanceFromWarehouse.end());
      }
      map<int, LSExpression>& vehicleEndMapTime =
          timesToWarehouses[vehicle.matrix_index()];
      if (vehicleEndMapTime.find(vehicle.end_index()) == vehicleEndMapTime.end()) {
        vector<int> timesToWarehouse;
        timesToWarehouse.reserve(problem.services_size());
        if (vehicle.end_index() > -1) {
          if (problem.matrices(vehicle.matrix_index()).time().size() > 0) {
            for (const auto& service : problem.services()) {
              timesToWarehouse.push_back(
                  problem.matrices(vehicle.matrix_index())
                      .time(service.matrix_index() * time_matrix_size +
                            vehicle.end_index()));
            }
          } else {
            for (const auto& service : problem.services()) {
              (void)service;
              timesToWarehouse.push_back(0);
            }
          }
        } else {
          timesToWarehouse.resize(problem.services_size());
          fill_n(timesToWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleEndMapTime[vehicle.end_index()] =
            model.array(timesToWarehouse.begin(), timesToWarehouse.end());
      }
      map<int, LSExpression>& vehicleEndMapTimeDistance =
          distanceToWarehouses[vehicle.matrix_index()];
      if (vehicleEndMapTimeDistance.find(vehicle.end_index()) ==
          vehicleEndMapTimeDistance.end()) {
        vector<int> distanceToWarehouse;
        distanceToWarehouse.reserve(problem.services_size());
        if (vehicle.end_index() > -1) {
          if (problem.matrices(vehicle.matrix_index()).distance().size() > 0) {
            for (const auto& service : problem.services()) {
              distanceToWarehouse.push_back(
                  problem.matrices(vehicle.matrix_index())
                      .distance(service.matrix_index() * distance_matrix_size +
                                vehicle.end_index()));
            }
          } else {
            for (const auto& service : problem.services()) {
              (void)service;
              distanceToWarehouse.push_back(0);
            }
          }
        } else {
          distanceToWarehouse.resize(problem.services_size());
          fill_n(distanceToWarehouse.begin(), problem.services_size(), 0);
        }
        vehicleEndMapTimeDistance[vehicle.end_index()] =
            model.array(distanceToWarehouse.begin(), distanceToWarehouse.end());
      }
      // cout << vehicle.id() << " " <<
      // vehicleStartMapTime[vehicle.start_index()].toString()
      //      << endl;
      // cout << vehicle.id() << " " << vehicleEndMapTime[vehicle.end_index()].toString()
      //      << endl;

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

      vehicle_ids_map_[(string)vehicle.id()] = k;
      k++;
    }
    vehicleStartIndicies =
        model.array(vehicleStartIndiciesVec.begin(), vehicleStartIndiciesVec.end());
    vehicleEndIndicies =
        model.array(vehicleEndIndiciesVec.begin(), vehicleEndIndiciesVec.end());

    vector<int> serviceTimeVec;
    vector<float> serviceLateMultiplierVec;
    vector<int> serviceSetUpDurationVec;
    vector<vector<float>> serviceQuantitiesVec;
    vector<float> serviceExclusionCostVec;
    vector<int> servicePriorityVec;
    vector<int> serviceMatrixIndexVec;
    vector<int> nbTWsVec;

    int s = 0;
    for (const auto& service : problem.services()) {
      serviceMatrixIndexVec.push_back(service.matrix_index());
      waitNextTWArray.addOperand(model.boolVar());
      serviceTimeVec.push_back(service.duration());
      serviceLateMultiplierVec.push_back(service.late_multiplier());
      serviceSetUpDurationVec.push_back(service.setup_duration());
      vector<float> serviceQuantityUnit;
      serviceQuantitiesVec.push_back(serviceQuantityUnit);
      for (int unit_i = 0; unit_i < service.quantities_size(); unit_i++) {
        serviceQuantitiesVec[s].push_back(service.quantities(unit_i));
      }
      serviceQuantitiesMatrix.addOperand(
          model.array(serviceQuantitiesVec[s].begin(), serviceQuantitiesVec[s].end()));
      if (service.exclusion_cost() > 0) {
        serviceExclusionCostVec.push_back(service.exclusion_cost() * CUSTOM_BIGNUM_COST);
      } else {
        serviceExclusionCostVec.push_back(1e6);
      }
      if (service.priority() != 4) {
        servicePriorityVec.push_back(service.priority());
      } else {
        servicePriorityVec.push_back(4);
      }
      vector<int> serviceTWStarts;
      vector<int> serviceTWEnds;
      vector<int> serviceTWAbsoluteEnds;

      if (service.time_windows_size() > 0) {
        for (const auto& tw : service.time_windows()) {
          serviceTWStarts.push_back(tw.start());
          serviceTWEnds.push_back(tw.end());
          if (static_cast<int>(tw.start()) > maxTwStarts) {
            maxTwStarts = tw.start();
          }
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
      serviceTwStartsArray.addOperand(
          model.array(serviceTWStarts.begin(), serviceTWStarts.end()));
      serviceTwStartsArray.setName("TW starts :");
      serviceTwEndsArray.addOperand(
          model.array(serviceTWEnds.begin(), serviceTWEnds.end()));
      twAbsoluteEndsArray.addOperand(
          model.array(serviceTWAbsoluteEnds.begin(), serviceTWAbsoluteEnds.end()));

      if (service.time_windows_size() == 0) {
        nbTWsVec.push_back(1);
      } else {
        nbTWsVec.push_back(service.time_windows_size());
      }
      service_ids_map_[(string)service.id()] = s;
      s++;
    }
    serviceMatrixIndex =
        model.array(serviceMatrixIndexVec.begin(), serviceMatrixIndexVec.end());
    serviceTime = model.array(serviceTimeVec.begin(), serviceTimeVec.end());
    serviceLateMultiplier =
        model.array(serviceLateMultiplierVec.begin(), serviceLateMultiplierVec.end());
    serviceSetUpDuration =
        model.array(serviceSetUpDurationVec.begin(), serviceSetUpDurationVec.end());
    serviceExclusionCost =
        model.array(serviceExclusionCostVec.begin(), serviceExclusionCostVec.end());
    servicePriority = model.array(servicePriorityVec.begin(), servicePriorityVec.end());
    nbTwsArray = model.array(nbTWsVec.begin(), nbTWsVec.end());
  }

  //   void precedenceConstraints(LSModel& model) {}

  void createModelAndSolve(localsolver_result::Result* result) {
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

    // all services must be satisfied by the vehicles
    LSExpression partition =
        model.partition(serviceSequences.begin(), serviceSequences.end());
    model.constraint(partition);

    LSExpression unassignedServices = serviceSequences[problem.vehicles_size()];
    LSExpression numberOfUnassignedServices = model.count(unassignedServices);
    vehiclesUsed[problem.vehicles_size()] = numberOfUnassignedServices > 0;

    vector<LSExpression> timeLeavingTheWarehouseConstraint(problem.vehicles_size());

    // Constraints for each vehicle
    int k = 0;
    for (const auto& vehicle : problem.vehicles()) {
      allVehicleIndices.push_back(k);
      LSExpression sequenceVehicle = serviceSequences[k];
      LSExpression c = model.count(sequenceVehicle);

      vehiclesUsed[k] = c > 0;
      vehiclesUsed[k].setName("vehicleUsed_" + to_string(k));

      routeFixedCost[k] = model.iif(
          vehiclesUsed[k], static_cast<lsint>(problem.vehicles(k).cost_fixed()), 0);

      lsint twStart = vehicle.has_time_window()
                          ? static_cast<lsint>(vehicle.time_window().start())
                          : 0;

      if (vehicle.shift_preference() == "force_start") {
        timeLeavingTheWarehouse[k] = model.intVar(twStart, twStart);
        startTimeVehicle[k] = timeLeavingTheWarehouse[k];
      } else if (maxTwStarts > 0) {
        timeLeavingTheWarehouse[k] = model.intVar(twStart, maxTwStarts);
        startTimeVehicle[k] = model.max(
            timeLeavingTheWarehouse[k],
            model.iif(c > 0,
                      model.max(twStart, model.at(twStartsArray, sequenceVehicle[0], 0) -
                                             timesFromWarehouses[vehicle.matrix_index()]
                                                                [vehicle.start_index()]
                                                                [sequenceVehicle[0]] -
                                             serviceSetUpDuration[sequenceVehicle[0]]),
                      0));
      } else {
        timeLeavingTheWarehouse[k] = model.intVar(twStart, twStart);
        startTimeVehicle[k] = timeLeavingTheWarehouse[k];
      }

      timeLeavingTheWarehouse[k].setName("timeLeavingTheWarehouse" + vehicle.id());

      timeLeavingTheWarehouseConstraint[k] =
          model.eq(timeLeavingTheWarehouse[k], twStart);
      model.constraint(timeLeavingTheWarehouseConstraint[k]);

      LSExpression distSelector = model.createLambdaFunction([&](LSExpression i) {
        return model.at(distanceMatrices[vehicle.matrix_index()], sequenceVehicle[i - 1],
                        sequenceVehicle[i]);
      });
      routeDistance[k] =
          model.sum(model.range(1, c), distSelector) +
          model.iif(c > 0,
                    distanceFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                          [sequenceVehicle[0]] +
                        distanceToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                            [sequenceVehicle[c - 1]],
                    0);
      routeDistanceCost[k] =
          (routeDistance[k] -
           model.iif(vehicle.free_approach(),
                     distanceFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                           [sequenceVehicle[0]],
                     0) -
           model.iif(vehicle.free_return(),
                     distanceToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                         [sequenceVehicle[c - 1]],
                     0)) *
          vehicle.cost_distance_multiplier();
      routeDistanceCost[k].setName("routeDistanceCost" + to_string(k));

      // End of each visit
      LSExpression endSelector;
      if (vehicle.rests_size() == 0) {
        endSelector = model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
          return nextStart(
                     sequenceVehicle[i],
                     model.iif(
                         i == 0,
                         startTimeVehicle[k] +
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
      } else {
        endSelector = model.createLambdaFunction([&](LSExpression i, LSExpression prev) {
          LSExpression tmp = travelAndSetupEndWithPause(
              vehicle, i, model.iif(i == 0, startTimeVehicle[k], prev), serviceSequences);
          return serviceEndWithPause(vehicle, i, tmp);
        });
      }
      endTime[k] = model.array(model.range(0, c), endSelector);

      LSExpression beginTimeSelector = model.createLambdaFunction([&](LSExpression i) {
        return endTime[k][i] - serviceTime[sequenceVehicle[i]];
      });
      beginTime[k] = model.array(model.range(0, c), beginTimeSelector);
      beginTime[k].setName("beginTime_" + problem.vehicles(k).id());

      LSExpression arrivalTimeSelector = model.createLambdaFunction([&](LSExpression i) {
        return beginTime[k][i] -
               model.iif(i < 1, serviceSetUpDuration[sequenceVehicle[i]],
                         model.iif(serviceMatrixIndex[sequenceVehicle[i]] ==
                                       serviceMatrixIndex[sequenceVehicle[i - 1]],
                                   0, serviceSetUpDuration[sequenceVehicle[i]]));
      });
      arrivalTime[k] = model.array(model.range(0, c), arrivalTimeSelector);
      arrivalTime[k].setName("arrivalTime" + problem.vehicles(k).id());

      LSExpression waitingTimeSelector = model.createLambdaFunction([&](LSExpression i) {
        return model.iif(
            i == 0,
            arrivalTime[k][i] -
                timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                   [sequenceVehicle[i]] -
                startTimeVehicle[k],
            model.max(0, arrivalTime[k][i] - endTime[k][i - 1] -
                             model.at(timeMatrices[vehicle.matrix_index()],
                                      sequenceVehicle[i - 1], sequenceVehicle[i])));
      });
      waitingTime[k] = model.sum(model.range(0, c), waitingTimeSelector);

      LSExpression endOfVehicleTimeWindow(model.array());
      endOfVehicleTimeWindow.addOperand(static_cast<lsint>(vehicle.time_window().end()));

      LSExpression pauseAfterLastServiceSelector =
          model.createLambdaFunction([&](LSExpression rest_index) {
            return model.iif(needsPause(Rest[k][rest_index], endTime[k][c - 1],
                                        endOfVehicleTimeWindow[k]),
                restDuration[k][rest_index], 0);
          });

      LSExpression pauseAfterLastService =
          model.sum(model.range(0, vehicle.rests_size()), pauseAfterLastServiceSelector);

      routeDuration[k] =
          model.iif(c > 0,
                    endTime[k][c - 1] + pauseAfterLastService +
                        timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                         [sequenceVehicle[c - 1]] -
                        startTimeVehicle[k],
                    0);

      routeDurationCost[k] = model.iif(
          c > 0,
          ((routeDuration[k] -
            model.iif(vehicle.free_approach(),
                      timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                                         [sequenceVehicle[0]],
                      0) -
            model.iif(vehicle.free_return(),
                      timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                       [sequenceVehicle[c - 1]],
                      0))) *
              vehicle.cost_time_multiplier(),
          0);

      routeDurationCost[k].setName("routeDurationCost" + to_string(k));

      timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()].setName(
          "timesToWarehouses of vehicle" + to_string(k));
      timesFromWarehouses[vehicle.matrix_index()][vehicle.start_index()].setName(
          "timesFromWarehouses of vehicle" + to_string(k));
      if (vehicle.duration() > 0 &&
          (!vehicle.has_time_window() ||
           vehicle.time_window().end() - vehicle.time_window().start() >
               vehicle.duration())) {
        LSExpression duration_constraint =
            routeDuration[k] <= static_cast<lsint>(vehicle.duration());
        duration_constraint.setName("duration_constraint_" + to_string(k));
        model.constraint(duration_constraint);
      }

      if (vehicle.has_time_window() && vehicle.time_window().end() < CUSTOM_MAX_INT) {
        if (vehicle.cost_late_multiplier() > 0) {
          LSExpression twEndsVehicleConstraintWithCostLateMultiplier = model.iif(
              c > 0,
              endTime[k][c - 1] +
                      timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                       [sequenceVehicle[c - 1]] <=
                  static_cast<lsint>(vehicle.time_window().end() +
                                     vehicle.time_window().maximum_lateness()),
              true);
          twEndsVehicleConstraintWithCostLateMultiplier.setName(
              "vehicle " + to_string(k) +
              "ending time window constraints with cost late multiplier");
          model.constraint(twEndsVehicleConstraintWithCostLateMultiplier);
        } else {
          LSExpression twEndsVehicleConstraint =
              model.iif(c > 0,
                        timesToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                                         [sequenceVehicle[c - 1]] +
                                endTime[k][c - 1] <=
                            static_cast<lsint>(vehicle.time_window().end()),
                        true);
          model.constraint(twEndsVehicleConstraint);
          twEndsVehicleConstraint.setName("vehicle " + to_string(k) +
                                          " ending time window constraints");
        }
      }
      LSExpression excessLatenessSelector =
          model.createLambdaFunction([&](LSExpression i) {
            // excess latenessCost is only calculated wrt the very last absolute end
            // because nextStart() can't return an infeasable intermediate time.
            return model.max(0, endTime[k][i] - serviceTime[sequenceVehicle[i]] -
                                    (model.at(twAbsoluteEndsArray, sequenceVehicle[i],
                                              nbTwsArray[sequenceVehicle[i]] - 1)));
          });
      excessLateness[k] = model.sum(model.range(0, c), excessLatenessSelector);

      LSExpression latenessSelector = model.createLambdaFunction([&](LSExpression i) {
        return serviceLateMultiplier[sequenceVehicle[i]] *
               model.max(0, beginTime[k][i] -
                                twEndSelect(sequenceVehicle[i], beginTime[k][i]));
      });

      // latenessTW = model.range(0, nbTwsArray[sequenceVehicle[i]]);
      latenessOfServicesOfVehicle[k] = model.array(model.range(0, c), latenessSelector);

      latenessCost[k] =
          model.iif(c > 0, model.sum(model.range(0, c), latenessSelector), 0);

      capacityConstraintsOfVehicle(k, sequenceVehicle, c);

      RelationBuilder(sequenceVehicle, unassignedServices);
      k++;
    }

    int s = 0;
    for (auto const& service : problem.services()) {
      vector<int64> incompatibleVehicleIndices;
      set_difference(allVehicleIndices.begin(), allVehicleIndices.end(),
                     service.vehicle_indices().begin(), service.vehicle_indices().end(),
                     back_inserter(incompatibleVehicleIndices));
      for (int64 incompatible_vehicle : incompatibleVehicleIndices) {
        model.constraint(!model.contains(serviceSequences[incompatible_vehicle], s));
      }
      s++;
    }

    LSExpression totalExcessLateness =
        model.sum(excessLateness.begin(), excessLateness.end());
    model.constraint(totalExcessLateness == 0);

    LSExpression totalLatenessCost =
        model.max(0, model.sum(latenessCost.begin(), latenessCost.end()));
    totalLatenessCost.setName("total Lateness Cost");

    LSExpression totalWaitingTime = model.sum(waitingTime.begin(), waitingTime.end());

    LSExpression exclusionCostCumulator = model.createLambdaFunction(
        [&](LSExpression i) { return serviceExclusionCost[unassignedServices[i]]; });

    LSExpression totalExclusionCost =
        model.sum(model.range(0, numberOfUnassignedServices), exclusionCostCumulator);
    totalExclusionCost.setName("total Exclusion Cost");

    // Total waitingTimebeforeLeavingWarehouse
    LSExpression totalWaitingTimeBeforeStart;
    totalWaitingTimeBeforeStart =
        model.sum(timeLeavingTheWarehouse.begin(), timeLeavingTheWarehouse.end());

    // Total vehicles used :
    nbVehiclesUsed = model.sum(vehiclesUsed.begin(), vehiclesUsed.end());

    // Total duration traveled :
    LSExpression totalDuration;
    totalDuration = model.sum(routeDuration.begin(), routeDuration.end());
    totalDuration.setName("totalDuration");

    // Total duration cost :
    LSExpression totalDurationCost;
    totalDurationCost = model.sum(routeDurationCost.begin(), routeDurationCost.end());
    totalDurationCost.setName("totalDurationCost");

    // Total distance traveled :
    LSExpression totalDistance;
    totalDistance = model.sum(routeDistance.begin(), routeDistance.end());
    totalDistance.setName("totalDistance");

    // Total distance cost :
    LSExpression totalDistanceCost;
    totalDistanceCost = model.sum(routeDistanceCost.begin(), routeDistanceCost.end());
    totalDistanceCost.setName("totalDistanceCost");

    // Total Fixed Costs
    LSExpression totalFixedCost;
    totalFixedCost = model.sum(routeFixedCost.begin(), routeFixedCost.end());
    totalFixedCost.setName("totalFixedCost");

    model.minimize(totalExclusionCost + totalDurationCost + totalDistanceCost +
                   totalLatenessCost + totalFixedCost);

    setInitialSolution(serviceSequences);

    firstAndSecondSolving(timeLeavingTheWarehouseConstraint);

    ParseSolution(result, serviceSequences, vehiclesUsed);

    // cout << " ----------------------Waiting times--------------------------  " << endl;
    // for (int v = 0; v < problem.vehicles_size(); v++) {
    //   cout << "waiting time service of " << problem.vehicles(v).id() << " "
    //        << waitingTime[v].getValue() << endl;
    //   cout << "waiting time before start" << problem.vehicles(v).id() << " "
    //        << timeLeavingTheWarehouse[v].getValue() << endl;
    // }

    if (problem.services_size() < 20) {
      for (int service_index = 0; service_index < problem.services_size();
           service_index++) {
        cout << " late multiplier of service " + to_string(service_index)
             << serviceLateMultiplier.getArrayValue().toString() << endl;
      }
      cout << model.toString() << endl;
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
        cout << "distanceFromWarehouses of : " << vehicle.id() << " "
             << distanceFromWarehouses[vehicle.matrix_index()][vehicle.start_index()]
                    .getArrayValue()
                    .toString()
             << endl;
        cout << "distanceToWarehouses of : " << vehicle.id() << " "
             << distanceToWarehouses[vehicle.matrix_index()][vehicle.end_index()]
                    .getArrayValue()
                    .toString()
             << endl;

        cout << " Vehicle Time Window Start :  " << vehicle.time_window().start() << endl;
        cout << " Vehicle Time Window Ends  :  " << vehicle.time_window().end() << endl;

        cout << "Vehicles Rests : "
             << Rest[IdIndex(vehicle.id(), vehicle_ids_map_)].getArrayValue().toString()
             << endl;
        cout << "Vehicles Rests Duration : "
             << restDuration[IdIndex(vehicle.id(), vehicle_ids_map_)]
                    .getArrayValue()
                    .toString()
             << endl;
      }

      cout << "number of time windows " << nbTwsArray.getArrayValue().toString() << endl;

      cout << " Services Time Window Starts Array :  "
           << serviceTwStartsArray.getArrayValue().toString() << endl;

      cout << " Services Time Window Ends Array :  "
           << serviceTwEndsArray.getArrayValue().toString() << endl;
      cout << " Services Time Window Absolute End Array :  "
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
          cout << IndexId(servicesCollection[i], service_ids_map_) << " ";
        }
        cout << endl;
      }

      for (int v = 0; v < problem.vehicles_size(); v++) {
        if (vehiclesUsed[v].getValue() != 1)
          continue;
        cout << "assigned service(s) to " << problem.vehicles(v).id() << " : ";
        LSCollection servicesCollection = serviceSequences[v].getCollectionValue();
        for (lsint i = 0; i < servicesCollection.count(); i++) {
          cout << IndexId(servicesCollection[i], service_ids_map_) << " ";
        }
        cout << endl;
      }
      cout << " ----------------------Number of Vehicles Used "
              "--------------------------  "
           << endl;
      cout << "nbVehicle Used : " << nbVehiclesUsed.getValue() << endl;

      cout << " ----------------------Begin times--------------------------  " << endl;
      for (int v = 0; v < problem.vehicles_size(); v++) {
        cout << "begin times service of " << problem.vehicles(v).id() << " "
             << beginTime[v].getArrayValue().toString() << endl;
      }

      cout << " ----------------------Arrival times--------------------------  " << endl;
      for (int v = 0; v < problem.vehicles_size(); v++) {
        cout << "Arrival times service of " << problem.vehicles(v).id() << " "
             << arrivalTime[v].getArrayValue().toString() << endl;
      }

      cout << " ----------------------End times--------------------------  " << endl;
      for (int v = 0; v < problem.vehicles_size(); v++) {
        cout << "end times service of " << problem.vehicles(v).id() << " "
             << endTime[v].getArrayValue().toString() << endl;
        cout << "total route duration :  " << problem.vehicles(v).id() << " "
             << routeDuration[v].getValue() << endl;
      }
      cout << " ----------------------Waiting times--------------------------  " << endl;
      for (int v = 0; v < problem.vehicles_size(); v++) {
        cout << "waiting time service of " << problem.vehicles(v).id() << " "
             << waitingTime[v].getValue() << endl;
        cout << "waiting time before start" << problem.vehicles(v).id() << " "
             << timeLeavingTheWarehouse[v].getValue() << endl;
      }
      cout << endl;
      cout << " -------------------- LATENESS -----------------------------------"
           << endl;
      for (int v = 0; v < problem.vehicles_size(); v++) {
        cout << "lateness " << problem.vehicles(v).id() << " "
             << latenessOfServicesOfVehicle[v].getArrayValue().toString() << endl;
      }

      cout << " total Absolute Latenness " << totalExcessLateness.getValue() << endl;
      cout << " total Latenness " << totalLatenessCost.getDoubleValue() << endl;
      cout << endl;
      for (int v_index = 0; v_index < problem.vehicles_size(); v_index++) {
        cout << latenessOfServicesOfVehicle[v_index].getArrayValue().toString() << endl;
      }

      cout << " -------------------- OBJECTIVE -----------------------------------"
           << endl;
      for (int i = 0; i < model.getNbObjectives(); i++) {
        cout << model.getObjective(i).toString() << endl;
      }
      cout << endl;
      cout << " ----------------------- END "
              "-------------------------------------------"
           << endl;
    }
  }
};

void readData(localsolver_vrp::Problem& problem) {
  const string filename = absl::GetFlag(FLAGS_instance_file);
  fstream input(filename, ios::in | ios::binary);
  if (!problem.ParseFromIstream(&input)) {
    cout << "Failed to parse protobuf." << endl;
  }

  for (const auto& vehicle : problem.vehicles()) {
    // cout << "cost waiting time multiplier : " <<
    // vehicle.cost_waiting_time_multiplier()
    //      << endl;
    // cout << "cost time multiplier : " << vehicle.cost_time_multiplier() << endl;
    if (vehicle.max_ride_distance()) {
      throw std::invalid_argument(" ERROR ======================= "
                                  " Max ride distance is not implemented yet");
    }
    if (vehicle.max_ride_time())
      throw std::invalid_argument(" ERROR ======================== "
                                  "Max ride time is not implemented yet");
    if (vehicle.value_matrix_index()) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "Value matrix is not implemented yet");
    }
    if (!equals(vehicle.cost_value_multiplier(), 0)) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "cost_value_multiplier is not implemented yet");
    }
    // if (vehicle.rests_size() > 0) {
    //   throw std::invalid_argument(" ERROR ======================= "
    //                               "rests are not implemented yet");
    // }
    if ((!vehicle.shift_preference().empty()) &&
        (vehicle.shift_preference() != "minimize_span" &&
         vehicle.shift_preference() != "force_start")) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "shift_preference is not implemented yet");
    }
    if (!equals(vehicle.additional_service(), 0)) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "additional service is not implemented yet");
    }
    if (!equals(vehicle.additional_setup(), 0)) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "additional setup is not implemented yet");
    }

    if (!equals(vehicle.coef_service(), 1)) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "coef_service is not implemented yet");
    }
    if (!equals(vehicle.coef_setup(), 1)) {
      throw std::invalid_argument(" ERROR ======================= "
                                  "coef_setup is not implemented yet");
    }
    // for (const auto& capacity : vehicle.capacities()) {
    //   if (capacity.overload_multiplier() > 0) {
    //     throw std::invalid_argument(" ERROR ======================= "
    //                                 "overload_multiplier is not implemented yet");
    //   }
    // }
  }

  for (const auto& service : problem.services()) {
    for (int quantity_index = 0; quantity_index < service.refill_quantities_size();
         quantity_index++) {
      if (service.refill_quantities(quantity_index)) {
        throw std::invalid_argument(" ERROR ======================= "
                                    " refill quantities are not implemented yet");
      }
    }
  }
};

int main(int argc, char** argv) {
  signal(SIGSEGV, handler);
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  localsolver_vrp::Problem problem;
  readData(problem);
  if (problem.services_size() < 10) {
    for (auto vehicle : problem.vehicles()) {
      vehicle.PrintDebugString();
    }

    for (auto service : problem.services()) {
      cout << "duration of service " << service.id() << " : " << service.duration()
           << endl;
      cout << "exclusion cost :" << service.exclusion_cost() << endl;
    }
  }
  localsolver_result::Result* result = new localsolver_result::Result;
  localsolver_VRP model(problem);

  model.createModelAndSolve(result);
  if (result != nullptr) {
    std::ofstream output(absl::GetFlag(FLAGS_solution_file),
                         std::ios::trunc | std::ios::binary);
    if (!result->SerializeToOstream(&output)) {
      cout << "Failed to write result.";
    }
    output.close();
  } else {
    std::cout << "No solution found..." << std::endl;
  }
  delete result;

  gflags::ShutDownCommandLineFlags();
  return 0;
}