#include <iostream>
#include <fstream>
#include "function.hpp"
#include "jacobian.hpp"
#include "LUDecompSolver.hpp"

int main(int argc, char **argv)
{   
    /*
    double x10, x20, v1, v2, k1, k2, m1, m2, l1, l2;
    if (argc < 11)
    {
        fprintf(
            stderr,
            "Usage: prog <x1> <x2> <v1> <v2>  <k1> <k2> <m1> <m2> <l1> <l2>\n"
        );
        return -1;
    }
    if (sscanf(argv[1], "%lf", &x10) < 1 ||
        sscanf(argv[2], "%lf", &x20) < 1 ||
        sscanf(argv[3], "%lf", &v1) < 1 ||
         sscanf(argv[4], "%lf", &v2) < 1 ||
        sscanf(argv[5], "%lf", &k1) < 1 ||
         sscanf(argv[6], "%lf", &k2) < 1 ||
        sscanf(argv[7], "%lf", &m1) < 1 ||
         sscanf(argv[8], "%lf", &m2) < 1 ||
        sscanf(argv[9], "%lf", &l1) < 1 || 
        sscanf(argv[10], "%lf", &l2) < 1) 
    {
      fprintf(stderr,
       "Usage: prog <x1> <x2> <v1> <v2>  <k1> <k2> <m1> <m2> <l1> <l2>\n");
      return -1;
    }
    Params p = {x10, x20, v1, v2, k1, k2, m1, m2, l1, l2};
    */
    // target param values to run out model to generate
    // data points to simulate measurement data 

    Params targetParams = {1, 1.4, 1, 1, 1, 1, 1, 1, 1, 1};
    // initial param values to start optimization
    Params optParams = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    FourColumnTable stepsData;
    double odeIntegrationStep = 0.001;
    size_t maxSteps = 60000;
    double noiseLevel = 0.2;

    //jacobian numberical calculation step
    double paramStep = 0.001;
    size_t gaussNewtonIterations = 100;
    // DOPRI8(array, p, steps, h, stepsData);
    // Euler(array, p, steps, h, stepsData);
    // total number of data sample is about maxSteps / sampleDistance
    size_t sampleDistance = 500;
    std::vector<int> solutionSteps; // = {0, 2000, 3000};
    
    for (int i = 0; i < maxSteps; i++) {
        if (i % sampleDistance == 0) {
            solutionSteps.push_back(i);
        }
    }
    
    FourColumnTable solutionOutput;
    DOPRI8(
      solutionOutput,
      targetParams, solutionSteps, odeIntegrationStep, stepsData);
    // sample some points from the stepsData
    FourColumnTable dataPoints;
    dataPoints.resize(solutionOutput.size());
    srand(time(NULL));
    for(auto solId = 0; solId < solutionOutput.size(); solId++)
    {
        for (auto i = 0; i<4; ++i){
            double noisy = 
                solutionOutput[solId][i] 
                + (rand() / double(RAND_MAX) - 0.5) * noiseLevel;
            dataPoints[solId][i] = noisy;
        }
    }

    writeSolutionToFiles(
        "data_points.txt", dataPoints,
        "data.txt", solutionSteps, stepsData);

    // here we use gauss-newton method to find parameters

    std::cout << "target params: " << std::endl;
    for (const auto &el : targetParams.p) {
      std::cout << el << ' ';
    }
    std::cout << '\n';
    std::cout << "optParams: " << std::endl;
    for (const auto &el : optParams.p) {
      std::cout << el << ' ';
    }
    std::cout << '\n';

    Matrix target(dataPoints.size() * 4, 1);
    size_t rowIdx = 0;
    for (const auto &row : dataPoints) {
      for (const auto &el : row) {
        target.matrix[rowIdx][0] = el;
        ++rowIdx;
      }
    }

    double sumOfSquares = INFINITY;
    double prevSumOfSquares = INFINITY;
    for (size_t ic = 0; ic < gaussNewtonIterations; ++ic) {
      prevSumOfSquares = sumOfSquares;
      Matrix resid = residual(
        target, optParams, odeIntegrationStep, solutionSteps, sumOfSquares);
      std::cout << "iteration " << ic << std::endl;
      std::cout << "targetParams" << std::endl;
      for (size_t i = 0; i < PARAM_COUNT; ++i) {
        std::cout << targetParams.p[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "optParams" << std::endl;
      for (size_t i = 0; i < PARAM_COUNT; ++i) {
        std::cout << optParams.p[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "sum of squares: " << sumOfSquares 
                << std::endl;
      if ((sumOfSquares > prevSumOfSquares) || (std::isnan(sumOfSquares))) {
        std::cout << "sum of squares increased or nan" << std::endl;
        break;
      }    
      Matrix jac = jacobianOfResidual(
        optParams, paramStep, odeIntegrationStep, solutionSteps
     );
    
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

        optParams.p[row] += deltaParams.matrix[row][0];
      }
      
      std::cout << std::endl;
    
    }

    return 0;
}
