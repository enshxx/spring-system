#ifndef JACOBIAN_H
#define JACOBIAN_H
#include "function.hpp"
#include "matrix.hpp"
#include "DOPRINew.hpp"
#include <sstream>

// #define WRITE_SOLUTIONS_TO_FILES

Matrix residual(Matrix const &target,
                Params const &p,
                double integrationStep,
                std::vector<int> const &timeSteps,
                double& sumOfSquares) {

  Matrix result(target.rows(), target.cols());

  sumOfSquares = 0.0;

  FourColumnTable stepsDebugData;
  FourColumnTable solutionPoints;
  DOPRI8(solutionPoints, p, timeSteps, integrationStep, stepsDebugData);

  for (auto solId = 0; solId < solutionPoints.size(); ++solId) {
    for (auto i = 0; i < 4; ++i) {
      auto rowId = 4 * solId + i;
      double r = target.matrix[rowId][0] - solutionPoints[solId][i];
      sumOfSquares += r * r;
      result.matrix[rowId][0] = r;
    }
  }
  return result;
}

Matrix jacobianOfResidual(Params const &p, double paramStep,
                          double integrationStep,
                          std::vector<int> const &timeSteps) {
  size_t DataPointCount = timeSteps.size();
  Matrix result(4 * DataPointCount, PARAM_COUNT);

  FourColumnTable stepsDebugData;
  for (int paramIdx = 0; paramIdx < PARAM_COUNT; ++paramIdx) {

    FourColumnTable statePlus, stateMinus;
    DOPRI8(statePlus, p.step(paramStep, static_cast<ParamIndices>(paramIdx)),
           timeSteps, integrationStep, stepsDebugData);
#ifdef WRITE_SOLUTIONS_TO_FILES
    std::stringstream ss;
    ss << "plus_data" << paramIdx  << ".txt";
    writeSolutionToFiles(
        "", statePlus,
        ss.str(), timeSteps, stepsDebugData);
#endif
    DOPRI8(stateMinus, p.step(-paramStep, static_cast<ParamIndices>(paramIdx)),
           timeSteps, integrationStep, stepsDebugData);
#ifdef WRITE_SOLUTIONS_TO_FILES
    ss.str("");
    ss << "minus_data" << paramIdx  << ".txt";
    writeSolutionToFiles(
        "", stateMinus,
        ss.str(), timeSteps, stepsDebugData);
#endif
    for (auto solId = 0; solId < statePlus.size(); ++solId) {
      for (auto j = 0; j < 4; ++j) {
        auto rowId = 4 * solId + j;
        result.matrix[rowId][paramIdx] =
            -1.0 * (statePlus[solId][j] - stateMinus[solId][j]) /
            (2 * paramStep);
      }
    }
  }
  return result;
}


//calculate Jacbian based on Euler method finite difference
// approximation for our ODE system
// params here assumed to be 
Matrix jacobianOfResidualEuler(Params const &p, double paramStep,
                          double integrationStep,
                          std::vector<int> const &timeSteps) {
  size_t DataPointCount = timeSteps.size();
  uint8_t velODEcount = 2;
  Matrix result(velODEcount * DataPointCount, PARAM_COUNT);
  FourColumnTable solution;
  FourColumnTable stepsDebugData;
  DOPRI8(solution, p, timeSteps, integrationStep, stepsDebugData);

  size_t velocityEqCount = 2;
  for (auto solId = 0; solId < solution.size(); ++solId) {
    // assume dx1/dparam_i, dx2/dparam_i is zero 
    // because they depend on dt*dt
    // so jacobian only calculated for 
    for (auto j = 0; j < velocityEqCount; ++j) {
      auto rowId = velocityEqCount * solId + j;
      double dt = integrationStep;
      double x1n = solution[solId][X1];
      double x2n = solution[solId][X2];
      result.matrix[rowId][K1] = -dt*(p[L1] - x1n)/p[M1];
      result.matrix[rowId][K2] = 
        -dt*(-p[L2] - x1n + x2n)/p[M1];
      result.matrix[rowId][M1] = 
        dt*(-p[K1]*(-p[L1] + x1n) + p[K2]*(-p[L2] - x1n + x2n))/(p[M1]*p[M1]);
      result.matrix[rowId][M2] = 0.0;
      result.matrix[rowId][L1] = -dt*p[K1]/p[M1];
      result.matrix[rowId][L2] = dt*p[K2]/p[M1];
    }
  }

  return result;
}


#endif