#include <iostream>
#include <fstream>
#include "jacobian.hpp"
#include "LUDecompSolver.hpp"
#include "jacobianPlusODE.hpp"


void printParams(std::__1::vector<double> &targetParams, std::__1::vector<double> &optParams)
{
  std::cout << "target params: " << std::endl;
  for (const auto &el : targetParams)
  {
    std::cout << el << ' ';
  }
  std::cout << '\n';
  std::cout << "optParams: " << std::endl;
  for (const auto &el : optParams)
  {
    std::cout << el << ' ';
  }
  std::cout << '\n';
}

int main(int argc, char **argv)
{   
// params legend  {v10 x10 v20 x20 K1 K2 M1 M2 L1 L2}
    std::vector<double> targetParams = 
        {1, 1, 1, 1.4, 1, 1, 1, 1, 1, 1};
    // initial param values to start optimization
    std::vector<double> optParams = {1, 0.7, 0.0, 1.3, 1, 1, 1, 1, 1, 1};
    
    double odeIntegrationStep = 0.001;
    size_t maxSteps = 60000;
    double noiseLevel = 0.2;
  
    size_t gaussNewtonIterations = 100;
    // DOPRI8(array, p, steps, h, stepsData);
    // Euler(array, p, steps, h, stepsData);
    // total number of data samples is about maxSteps / sampleDistance
    size_t sampleSteps = 600;
    
    // sample some points from the stepsData
    Matrix dStatedTheta;
    Matrix stateDataPoints;
    Matrix generatedMeasurements(sampleSteps*STATE_SIZE, 1);
    solveStatePlusJacODEs(
      targetParams,
      maxSteps,
      odeIntegrationStep,
      sampleSteps,
      stateDataPoints,
      dStatedTheta
    );
    printMatrixToFile(dStatedTheta, "jacobian.txt");
    srand(time(NULL));
    for(auto solId = 0; solId < stateDataPoints.rows(); solId++)
    {
        for (auto i = 0; i<STATE_SIZE; ++i){
            double noisy = 
                stateDataPoints.matrix[solId][0] 
                + (rand() / double(RAND_MAX) - 0.5) * noiseLevel;
            generatedMeasurements.matrix[solId][0] = noisy;
        }
    }
    printMatrixToFile(generatedMeasurements, "data_points.txt");


    // here we use gauss-newton method to find parameters

    printParams(targetParams, optParams);

    double sumOfSquares = INFINITY;
    double prevSumOfSquares = INFINITY;
  
    for (size_t ic = 0; ic < gaussNewtonIterations; ++ic) {
      prevSumOfSquares = sumOfSquares;
      Matrix resid(generatedMeasurements.rows(), 1);
      Matrix stateSolution;
      Matrix jac;
      solveStatePlusJacODEs(optParams, maxSteps, odeIntegrationStep, 
        sampleSteps, stateSolution, jac);
      // multipy jacobian by -1 
      multiplyMatrixByNumber(jac, -1.0);
      subtractMatrices(generatedMeasurements, stateSolution, resid);
      // printMatrixToFile(resid, "resid.txt");
      Matrix residTransposed = transpose(resid);
      Matrix sumOfSquaresMatrix;
      multiplyMatrices(residTransposed, resid, sumOfSquaresMatrix);
      sumOfSquares = sumOfSquaresMatrix.matrix[0][0];
      std::cout << "iteration " << ic << std::endl;
      printParams(targetParams, optParams);
      std::cout << std::endl;
      std::cout << "sum of squares: " << sumOfSquares
                << std::endl;
      if (/*(sumOfSquares > 2*prevSumOfSquares)* || */ (std::isnan(sumOfSquares))) {
        std::cout << "sum of squares is nan, exitting" << std::endl;
        break;
      }    
      // printMatrixToFile(jac, "jacobian.txt");
      Matrix jacTransposed = transpose(jac);
      Matrix jacTMulByRes;
      multiplyMatrices(jacTransposed, resid, jacTMulByRes);

      Matrix jacTByJac;
      multiplyMatrices( jacTransposed, jac, jacTByJac);
      multiplyMatrixByNumber(jacTMulByRes, -1.0);
      // solve linear system jacTByJac * x = jacTMulByRes
      int n = jacTMulByRes.rows();
      Matrix L, U;
      luDecomposition(jacTByJac, n, L, U);
      Matrix deltaParams(n, 1);
      solveLinearSystem(L, U, jacTMulByRes, deltaParams, n);

      for (size_t row = 0; row < deltaParams.rows(); ++row) {

        optParams[row] += deltaParams.matrix[row][0];
      }    
    }
  
    return 0;
}
