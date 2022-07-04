#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <lsarray.h>
#include <lscollection.h>
#include <lsexpression.h>
#include <lssolution.h>
#include <ostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "./localsolver_result.pb.h"
#include "./localsolver_vrp.pb.h"
#include "localsolver.h"
#include "ortools/base/commandlineflags.h"

bool compareFloatNum(double a, double b) {
  // Correct method to compare
  // floating-point numbers
  if (abs(a - b) < 1e-9) {
    return true;
  } else {
    return false;
  }
}