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

#endif