#include "jacobianPlusODE.hpp"
#include "matrix.hpp"


int main() {

    std::vector<double> params = 
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double odeIntegrationStep = 0.001;
    size_t maxSteps = 60000;
    Matrix state, diffStateByParams;
    solveStatePlusJacODEs(
        params,
        maxSteps, 
        odeIntegrationStep,
        600,
        state,
        diffStateByParams); 
    printMatrixToFile(state, "State.txt");
    printMatrixToFile(diffStateByParams, "diffStateByParams.txt");
    
    return 0;  
}