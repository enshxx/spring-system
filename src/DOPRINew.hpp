#ifndef DOPRINEW_H
#define DOPRINEW_H

#include <cmath>
#include <iostream>
#include <vector>

#include "function.hpp"
#include "jacobianPlusODE.hpp"
#include "odesolvers.hpp"


// Define the system of 4 ODEs
void fnew(
  double t, double y[], double dydt[],
   const double * p, MatrixData& mData) {
  Params pp;
  for (auto i = 0; i < PARAM_COUNT; i++) {
    pp.p[i] = p[i];
  }
  StateVector derivatives = f({y[0], y[1], y[2], y[3]}, pp);
  for (int i = 0; i < 4; i++) {
      dydt[i] = derivatives[i];
  }
}

void DOPRI8(
    FourColumnTable &array,
    Params const& p,
    std::vector<int> const &steps,
    double h,
    std::vector<StateVector > &stepPoints) { //DOPNEW
  const int n = 4; // number of ODEs
  double y0[n] = {p[X10], p[X20], p[V10], p[V20]}; // initial conditions
  double t0 = 0.0; // initial time
 
  double t = 0;
  double y[n];
  for (int i = 0; i < n; i++) {
      y[i] = y0[i];
  }
  if (steps.empty()) { return; }
  array.resize(steps.size());
  size_t maxSteps = steps[steps.size()-1];
  stepPoints.resize(maxSteps);

  size_t currentTimeStep = 0;
  MatrixData m;
  for (size_t i = 0; i <= maxSteps; i++) {
    if (i == steps[currentTimeStep]) {
      array[currentTimeStep++] = {y[0], y[1], y[2], y[3]};
    } 
    
    dormandPrince(t, y0, h, t, y, p.p.data(), fnew, 4, m);
    stepPoints[i] = {y[0], y[1], y[2], y[3]};
    for (int i = 0; i < n; i++) {
      y0[i] = y[i];
    }
    
  }

}

#endif