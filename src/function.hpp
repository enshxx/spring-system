#ifndef FUNCTION_H
#define FUNCTION_H
#include <array>
#include <vector>
#include <string>
#include <fstream>


const int X1 = 0;
const int X2 = 1;
const int V1 = 2;
const int V2 = 3;

#define PARAM_COUNT 10

enum ParamIndices {
    X10 = 0,
    X20 = 1,
    V10 = 2,
    V20 = 3,
    K1 = 4,
    K2 = 5,
    M1 = 6,
    M2 = 7,
    L1 = 8,
    L2 = 9
};

using FourColumnTable = std::vector<std::array<double, 4> >;
using StateVector = std::array<double, 4>;

struct Params{
    std::array<double, 10> p;
    Params step(double h, ParamIndices index) const {
        Params result = {p};  
        result.p[index] += h;
        return result;
    }
    double operator[](ParamIndices index) const { return p[index]; }
};

std::array<double, 4> f(
    std::array<double, 4> state,
    Params const& p){
    std::array<double, 4> derivative;
    double k1 = p[K1];
    double k2 = p[K2];
    double m1 = p[M1];
    double m2 = p[M2];
    double l1 = p[L1];
    double l2 = p[L2];
    derivative[X1] = state[V1];
    derivative[X2] = state[V2];
    derivative[V1] = ((-k1 * (state[X1] - l1)) + (k2 * (state[X2] - state[X1] - l2))) / m1;
    derivative[V2] = (-k2 * (state[X2] - state[X1] - l2)) / m2;
    return derivative;
}

void Euler(
    std::array<double, 4> &array,
    Params const& p,
    int steps,
    double h,
    std::vector<std::array<double, 4> > &stepPoints) {

    stepPoints.resize(steps);
    for (int i = 0; i < steps; i++) {
        std::array<double, 4> derivative = f(array, p);
        array[X1] += h * derivative[X1];
        array[X2] += h * derivative[X2];
        array[V1] += h * derivative[V1];
        array[V2] += h * derivative[V2];
        stepPoints[i] = array;
    }
}

void writeSolutionToFiles(std::string const &pointsName,
                          FourColumnTable &dataPoints,
                          std::string const &solName,
                          std::vector<int> const &solutionSteps,
                          FourColumnTable &stepsData) {
  if (pointsName.length() > 0) {
    std::ofstream outDataPoints(pointsName);
    if (!outDataPoints.is_open()) {
      fprintf(stderr, "Can't open file %s\n", pointsName.c_str());
      return;
    }
    size_t timeStepId = 0;
    for (const auto &row : dataPoints) {
      for (const auto &el : row) {
        outDataPoints << el << ' ';
      }
      outDataPoints << solutionSteps[timeStepId++] << ' ';
      outDataPoints << std::endl;
    }
    outDataPoints.close();
  }
  std::ofstream out(solName);
  if (!out.is_open()) {
    fprintf(stderr, "Can't open file %s\n", solName.c_str());
    return;
  }
  for (const auto &row : stepsData) {
    for (const auto &el : row) {
      out << el << ' ';
    }
    out << '\n';
  }
  out.close();
}

#endif