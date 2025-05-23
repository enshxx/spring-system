#include <iostream>
#include <fstream>
#include <random>
#include "jacobian.hpp"
#include "linearsystemsolver.hpp"
#include "jacobianPlusODE.hpp"
#include "qrdecompsolver.hpp"


void printParams(std::vector<double> &targetParams, std::vector<double> &optParams)
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

double generateNoise(double sigma = 1.0)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> dist(0.0, sigma);
    return dist(gen);
}

int main(int argc, char **argv)
{   
// params legend  {v10 x10 v20 x20 M1 M2 L1 L2 K1 K2}
    std::vector<double> targetParams = 
        {-1, 2.5, 1, 1.4, 5, 3, 2, 5, 1, 3};
    // initial param values to start optimization
    std::vector<double> optParams = {-1, 2.4, 1.4, 1, 5.4, 3.4, 2, 5, 1, 3};
    
    double odeIntegrationStep = 0.001;
    double noiseLevel = 0.2;
  
    size_t gaussNewtonIterations = 100;
    // DOPRI8(array, p, steps, h, stepsData);
    // Euler(array, p, steps, h, stepsData);
    // total number of data samples is about maxSteps / sampleDistance
    size_t sampleSteps = 1000;
    // sample some points from the stepsData
    Matrix dStatedTheta;
    Matrix stateDataPoints;
    Matrix generatedMeasurements(sampleSteps*STATE_SIZE, 1);
    solveStatePlusJacODEs(
      targetParams,
      odeIntegrationStep,
      sampleSteps,
      stateDataPoints,
      dStatedTheta
    );
    printMatrixToFile(dStatedTheta, "jacobian.txt");
    // srand(time(NULL));
    for(auto solId = 0; solId < stateDataPoints.rows(); solId++)
    {
        for (auto i = 0; i<STATE_SIZE; ++i){
            double noisy = 
                stateDataPoints.matrix[solId][0]
                + generateNoise(noiseLevel);
            generatedMeasurements.matrix[solId][0] = noisy;
        }
    }
    printMatrixToFile(generatedMeasurements, "data_points.txt");


    // here we use gauss-newton method to find parameters

    printParams(targetParams, optParams);

    double sumOfSquares = INFINITY;
    double prevSumOfSquares = INFINITY;

    std::vector<double> prevOptParams = optParams;
    // vector regualarization values
    // values of each parameters are initially the same
    // but probaly can be adjusterd depending on how well 
    // given parameter is estimated
  
    for (size_t ic = 0; ic < gaussNewtonIterations; ++ic) {
      prevSumOfSquares = sumOfSquares;
      Matrix resid(generatedMeasurements.rows(), 1);
      Matrix stateSolution;
      Matrix jac;
      solveStatePlusJacODEs(
        optParams, odeIntegrationStep, 
        sampleSteps, stateSolution,
        jac);

      // printMatrixToFile(jac, "jacobian.txt");

      subtractMatrices(generatedMeasurements, stateSolution, resid);
      // printMatrixToFile(resid, "resid.txt");
      Matrix residTransposed(resid.cols(), resid.rows());
      transpose(resid, residTransposed);
      Matrix sumOfSquaresMatrix(
        residTransposed.rows(),resid.cols());
      multiplyMatrices(residTransposed, resid, sumOfSquaresMatrix);
      sumOfSquares = sumOfSquaresMatrix.matrix[0][0];
      std::cout << "iteration " << ic << std::endl;
      printParams(targetParams, optParams);
      std::cout << std::endl;
      std::cout << "sum of squares: " << sumOfSquares
                << std::endl;
      if (/*(sumOfSquares > 2*prevSumOfSquares)* || */ (std::isnan(sumOfSquares))) {
        std::cout << "sum of squares is nan, exitting" << std::endl;
        optParams = prevOptParams;
        break;
      }    
     
      Matrix jacTransposed(jac.cols(), jac.rows());
      transpose(jac, jacTransposed);

      Matrix jacTMulByRes(
        jacTransposed.rows(), resid.cols());
      multiplyMatrices(jacTransposed, resid, jacTMulByRes);

      Matrix jacTByJac(
        jacTransposed.rows(), jac.cols());
      multiplyMatrices(jacTransposed, jac, jacTByJac);

      // printMatrixToFile(jacTByJac, "jacTByJac.txt");
      
      // solve linear system jacTByJac * x = jacTMulByRes
      int n = jacTMulByRes.rows();
      Matrix deltaParams(n, 1);
      for (size_t row = 0; row < deltaParams.rows(); ++row) {
        deltaParams.matrix[row][0] = 0.0;
      }
      // solveLinearSystemLU( jacTByJac, jacTMulByRes, deltaParams);
      solveLinearSystemQR(jacTByJac, jacTMulByRes, deltaParams);
      prevOptParams = optParams;
      for (size_t row = 0; row < deltaParams.rows(); ++row) {

        optParams[row] += deltaParams.matrix[row][0];
      }    
    }

    solveStatePlusJacODEs(
      optParams,
      odeIntegrationStep,
      sampleSteps,
      stateDataPoints,
      dStatedTheta
    );
    printMatrixToFile(stateDataPoints, "opt_state_data_points.txt");
    // calculate solution for latest params

    return 0;
}