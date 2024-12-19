
#define STATE_SIZE 4
#define PARAM_SIZE 10

#include <array>
#include "matrix.hpp"
#include "DOPRINew.hpp"

// param names
constexpr size_t theta1 = 0;
constexpr size_t theta2 = 1;
constexpr size_t theta3 = 2;
constexpr size_t theta4 = 3;
constexpr size_t theta5 = 4;
constexpr size_t theta6 = 5;
constexpr size_t theta7 = 6;
constexpr size_t theta8 = 7;
constexpr size_t theta9 = 8;
constexpr size_t theta10 = 9;

constexpr size_t v10 = 0;
constexpr size_t x10 = 1;
constexpr size_t v20 = 2;
constexpr size_t x20 = 3;
constexpr size_t k1 = 4;
constexpr size_t k2 = 5;
constexpr size_t l1 = 6;
constexpr size_t l2 = 7;
constexpr size_t m1 = 8;
constexpr size_t m2 = 9;

// state names
constexpr size_t v1 = 0;
constexpr size_t x1 = 1;
constexpr size_t v2 = 2;
constexpr size_t x2 = 3;

constexpr size_t x_1 = 0;
constexpr size_t x_2 = 1;
constexpr size_t x_3 = 2;
constexpr size_t x_4 = 3;
    


void diffs(
    double t, 
    double* statePlusStateDiffTheta,
    double* diffs,
    const double* p
)
{
    // F - state diffs by t
    double *s = statePlusStateDiffTheta;
    double* stateDiffTheta = statePlusStateDiffTheta + STATE_SIZE;
    diffs[v1] = 
    (-p[theta5]*(-p[theta7] + s[x_2]) + p[theta6]*(-p[theta8] - s[x_2] + s[x_4]))/p[theta9];
    diffs[x1] = s[x_1];
    diffs[v2] = -p[theta6]*(-p[theta8] - s[x_2] + s[x_4])/p[theta10];
    diffs[x2] = s[x_3];
  
    Matrix FdiffState(STATE_SIZE, STATE_SIZE);
    FdiffState.matrix[0][0] = 0.0;
    FdiffState.matrix[0][1] = (-p[theta5] - p[theta6])/p[theta9];
    FdiffState.matrix[0][2] = 0.0;
    FdiffState.matrix[0][3] = p[theta6]/p[theta9];

    FdiffState.matrix[1][0] = 1;
    FdiffState.matrix[1][1] = 0;
    FdiffState.matrix[1][2] = 0;
    FdiffState.matrix[1][3] = 0;

    FdiffState.matrix[2][0] = 0;
    FdiffState.matrix[2][1] = p[theta6]/p[theta10];
    FdiffState.matrix[2][2] = 0;
    FdiffState.matrix[2][3] = -p[theta6]/p[theta10];

    FdiffState.matrix[3][0] = 0;
    FdiffState.matrix[3][1] = 0;
    FdiffState.matrix[3][2] = 1;
    FdiffState.matrix[3][3] = 0;

    Matrix FdiffParam(STATE_SIZE, PARAM_SIZE);
    FdiffParam.matrix[0][0] = 0.0;
    FdiffParam.matrix[0][1] = 0.0;    
    FdiffParam.matrix[0][2] = 0.0;
    FdiffParam.matrix[0][3] = 0.0;
    FdiffParam.matrix[0][4] = (p[theta7] - s[x_2])/p[theta9];
    FdiffParam.matrix[0][5] = (-p[theta8] - s[x_2] + s[x_4])/p[theta9];
    FdiffParam.matrix[0][6] = p[theta5]/p[theta9];
    FdiffParam.matrix[0][7] = -p[theta6]/p[theta9];
    FdiffParam.matrix[0][8] = 
        -(-p[theta5]*(-p[theta7] + s[x_2]) +
         p[theta6]*(-p[theta8] - s[x_2] + s[x_4]))/(p[theta9]*p[theta9]);
    FdiffParam.matrix[0][9] = 0.0;
       
    FdiffParam.matrix[1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    FdiffParam.matrix[2] = {
        0, 0, 0, 0, 0,
        -(-p[theta8] - s[x_2] + s[x_4])/p[theta10],
        0,
        p[theta6]/p[theta10],
        0,
        p[theta6]*(-p[theta8] - s[x_2] + s[x_4])/p[theta10]*p[theta10]
    };
    FdiffParam.matrix[3] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    Matrix stateDiffThetaDot(STATE_SIZE, PARAM_SIZE);

    Matrix tmp(STATE_SIZE, PARAM_SIZE);

    Matrix stateDiffThetaMat = Matrix(STATE_SIZE, PARAM_SIZE);
    for (auto row=0; row<STATE_SIZE; ++row) {
        stateDiffThetaMat.matrix[row].assign(
            stateDiffTheta + row*PARAM_SIZE, 
            stateDiffTheta + row*PARAM_SIZE+PARAM_SIZE
        );
    }
    multiplyMatrices(FdiffState, stateDiffThetaMat, tmp);
    addMatrices(tmp, FdiffParam, stateDiffThetaDot);

    for (size_t i = 0; i < STATE_SIZE; ++i) {
        for (size_t j = 0; j < PARAM_SIZE; ++j) {
            diffs[4 + i*PARAM_SIZE + j] = stateDiffThetaDot.matrix[i][j];
        }
    }

}


void solveStatePlusJacODEs( 
    std::vector<double> const & params,
    size_t solutionSteps, double odeIntegrationStep,
    size_t sampleSteps,
    Matrix &state,
    Matrix &diffStateByParams
) {
    
    state.matrix.resize(sampleSteps*STATE_SIZE);
    for (auto &s : state.matrix) {
        s.resize(1);
    }
    diffStateByParams.matrix.resize(sampleSteps*STATE_SIZE);
    for (auto &s : diffStateByParams.matrix) {
        s.resize(PARAM_SIZE);
    }

    // set initial values
    std::array<double, STATE_SIZE + STATE_SIZE*PARAM_SIZE> curSol = {
        // state at t = 0
        params[v10], params[x10], params[v20], params[x20],
        // diffStatebyTheta at t=0
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0
    };
    double t = 0;
    size_t sampleIncrement = solutionSteps / sampleSteps;
    size_t sampleId = 0;
    double debugDiffs[STATE_SIZE + STATE_SIZE*PARAM_SIZE];
    for (auto tStep = 0; tStep < solutionSteps; ++tStep) {
        
        dormandPrince(
            t,
            curSol.data(),
            odeIntegrationStep,
            t,
            curSol.data(),
            params.data(),
            diffs,
            STATE_SIZE + STATE_SIZE*PARAM_SIZE
        );
        
        if (tStep % sampleIncrement == 0) {

            double* p = curSol.data() + STATE_SIZE;
            for(auto i = 0; i < STATE_SIZE; ++i) {
                state.matrix[sampleId*STATE_SIZE+i][0] = curSol[i];
                diffStateByParams.matrix[sampleId*STATE_SIZE + i].assign(
                        p,  p+PARAM_SIZE
                    );
                p += PARAM_SIZE;
            }
            
            ++sampleId;
        }
    }
}

// Matrix residual();
