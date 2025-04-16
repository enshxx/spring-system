#ifndef ODESOLVERS_H
#define ODESOLVERS_H

#include "matrixdata.hpp"

void eulerSolver(double t0, double y0[], double h, double &t,
                 double y[], const double* p,
                 void (*dydt)(double t, double y[], double dydt[],
                  const double * p),
                 double* debugDiffs,
                 int n=4) {
    double k1[n];

    // Compute derivative at initial point
    dydt(t0, y0, k1, p);
    for (int i = 0; i < n; i++) {
        debugDiffs[i] = k1[i];
    }
    // Update each component
    for (int i = 0; i < n; i++) {
        y[i] = y0[i] + h * k1[i];
    }

    // Update time
    t = t0 + h;
}

// Define the Dormand-Prince method
void dormandPrince(double t0, double y0[], double h, double &t,
                  double y[],  const double*  p,
                  void (* dydt)(double t, double y[], double dydt[],
                    const double * p, MatrixData& mData),
                  int n,
                  MatrixData& mD
                   ) {
  
  double k1[n], k2[n], k3[n], k4[n], k5[n], k6[n];
  double yNew[n];

  // Compute k1
  dydt(t0, y0, k1, p, mD);

  // Compute k2
  double t2 = t0 + 1.0 / 5.0 * h;
  double y2[n];
  for (int i = 0; i < n; i++) {
    y2[i] = y0[i] + 1.0 / 5.0 * h * k1[i];
  }
  dydt(t2, y2, k2, p, mD);

  // Compute k3
  double t3 = t0 + 3.0 / 10.0 * h;
  double y3[n];
  for (int i = 0; i < n; i++) {
    y3[i] = y0[i] + 3.0 / 40.0 * h * k1[i] + 9.0 / 40.0 * h * k2[i];
  }
  dydt(t3, y3, k3, p, mD);

  // Compute k4
  double t4 = t0 + 4.0 / 5.0 * h;
  double y4[n];
  for (int i = 0; i < n; i++) {
    y4[i] = y0[i] + 44.0 / 45.0 * h * k1[i] - 56.0 / 15.0 * h * k2[i] +
            32.0 / 9.0 * h * k3[i];
  }
  dydt(t4, y4, k4, p, mD);

  // Compute k5
  double t5 = t0 + 8.0 / 9.0 * h;
  double y5[n];
  for (int i = 0; i < n; i++) {
    y5[i] = y0[i] + 19372.0 / 6561.0 * h * k1[i] -
            25360.0 / 2187.0 * h * k2[i] + 64448.0 / 6561.0 * h * k3[i] -
            212.0 / 729.0 * h * k4[i];
  }
  dydt(t5, y5, k5, p, mD);

  // Compute k6
  double t6 = t0 + h;
  double y6[n];
  for (int i = 0; i < n; i++) {
    y6[i] = y0[i] + 9017.0 / 3168.0 * h * k1[i] - 355.0 / 33.0 * h * k2[i] +
            46732.0 / 5247.0 * h * k3[i] + 49.0 / 176.0 * h * k4[i] -
            5103.0 / 18656.0 * h * k5[i];
  }
  dydt(t6, y6, k6, p, mD);

  // Compute yNew
  for (int i = 0; i < n; i++) {
    yNew[i] = y0[i] + 35.0 / 384.0 * h * k1[i] + 0.0 * h * k2[i] +
              500.0 / 1113.0 * h * k3[i] + 125.0 / 192.0 * h * k4[i] -
              2187.0 / 6784.0 * h * k5[i] + 11.0 / 84.0 * h * k6[i];
  }

  // Update t and y
  t = t0 + h;
  for (int i = 0; i < n; i++) {
    y[i] = yNew[i];
  }
}

#endif