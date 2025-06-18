#ifndef JACOBIANPLUSODE_H
#define JACOBIANPLUSODE_H


#include <array>
#include "matrix.hpp"
#include "odesolvers.hpp"
#include "matrixdata.hpp"

constexpr size_t dv1dt = 0;
constexpr size_t dx1dt = 1;
constexpr size_t dv2dt = 2;
constexpr size_t dx2dt = 3;
constexpr size_t dv1dv10 = 4;
constexpr size_t dv1dx10 = 5;
constexpr size_t dv1dv20 = 6;
constexpr size_t dv1dx20 = 7;
constexpr size_t dv1dm1 = 8;
constexpr size_t dv1dm2 = 9;
constexpr size_t dx1dv10 = 10;
constexpr size_t dx1dx10 = 11;
constexpr size_t dx1dv20 = 12;
constexpr size_t dx1dx20 = 13;
constexpr size_t dx1dm1 = 14;
constexpr size_t dx1dm2 = 15;
constexpr size_t dv2dv10 = 16;
constexpr size_t dv2dx10 = 17;
constexpr size_t dv2dv20 = 18;
constexpr size_t dv2dx20 = 19;
constexpr size_t dv2dm1 = 20;
constexpr size_t dv2dm2 = 21;
constexpr size_t dx2dv10 = 22;
constexpr size_t dx2dx10 = 23;
constexpr size_t dx2dv20 = 24;
constexpr size_t dx2dx20 = 25;
constexpr size_t dx2dm1 = 26;
constexpr size_t dx2dm2 = 27;

// constexpr size_t dv1dt = 0;
// constexpr size_t dx1dt = 1;
// constexpr size_t dv2dt = 2;
// constexpr size_t dx2dt = 3;
// constexpr size_t dv1dv10 = 4;
// constexpr size_t dv1dx10 = 5;
// constexpr size_t dv1dv20 = 6;
// constexpr size_t dv1dx20 = 7;
// constexpr size_t dv1dl1 = 8;
// constexpr size_t dv1dl2 = 9;
// constexpr size_t dv1dm1 = 10;
// constexpr size_t dv1dm2 = 11;
// constexpr size_t dx1dv10 = 12;
// constexpr size_t dx1dx10 = 13;
// constexpr size_t dx1dv20 = 14;
// constexpr size_t dx1dx20 = 15;
// constexpr size_t dx1dl1 = 16;
// constexpr size_t dx1dl2 = 17;
// constexpr size_t dx1dm1 = 18;
// constexpr size_t dx1dm2 = 19;
// constexpr size_t dv2dv10 = 20;
// constexpr size_t dv2dx10 = 21;
// constexpr size_t dv2dv20 = 22;
// constexpr size_t dv2dx20 = 23;
// constexpr size_t dv2dl1 = 24;
// constexpr size_t dv2dl2 = 25;
// constexpr size_t dv2dm1 = 26;
// constexpr size_t dv2dm2 = 27;
// constexpr size_t dx2dv10 = 28;
// constexpr size_t dx2dx10 = 29;
// constexpr size_t dx2dv20 = 30;
// constexpr size_t dx2dx20 = 31;
// constexpr size_t dx2dl1 = 32;
// constexpr size_t dx2dl2 = 33;
// constexpr size_t dx2dm1 = 34;
// constexpr size_t dx2dm2 = 35;

constexpr size_t v10 = 0;
constexpr size_t x10 = 1;
constexpr size_t v20 = 2;
constexpr size_t x20 = 3;
constexpr size_t m1 = 4;
constexpr size_t m2 = 5;
constexpr size_t l1 = 6;
constexpr size_t l2 = 7;
constexpr size_t k1 = 8;
constexpr size_t k2 = 9;
    

void calculateDiffs(
    double t, 
    double* stateDiff,
    double* diffs,
    const double* p
)
{
    
    double _dv1dx1 = -(p[k1] / p[m1] + p[k2] / p[m1]);
    double _dv1dx2 = p[k2] / p[m1];
    // double _dv1dk1 = -(stateDiff[dx1dt] - p[l1]) / p[m1];
    // double _dv1dk2 = (stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]) / p[m1];
    // double _dv1dl1 = p[k1] / p[m1];
    // double _dv1dl2 = -p[k2] / p[m1];
    double _dv1dm1 = p[k1] * (stateDiff[dx1dt] - p[l1]) / (p[m1] * p[m1]) - 
        p[k2] * (stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]) / (p[m1] * p[m1]);

    double _dv2dx1 = p[k2] / p[m2];
    double _dv2dx2 = -p[k2] / p[m2];
    // double _dv2dk2 = -(stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]) / p[m2];
    // double _dv2dl2 = p[k2] / p[m2];
    double _dv2dm2 = p[k2] * (stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]) / (p[m2] * p[m2]);

    for (auto i = 0; i < STATE_SIZE / 2;  ++i)
        for (auto j = 0; j < PARAM_SIZE; ++j)
            diffs[STATE_SIZE + 2 * PARAM_SIZE * i + PARAM_SIZE + j] = stateDiff[STATE_SIZE + 2 * PARAM_SIZE * i + j];

    diffs[dv1dt] = -p[k1] / p[m1] * (stateDiff[dx1dt] - p[l1]) + 
        p[k2] / p[m1] * (stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]);
    diffs[dx1dt] = stateDiff[dv1dt];
    diffs[dv2dt] = -p[k2] / p[m2] * (stateDiff[dx2dt] - stateDiff[dx1dt] - p[l2]);
    diffs[dx2dt] = stateDiff[dv2dt];

    diffs[dv1dv10] = _dv1dx1 * stateDiff[dx1dv10] + _dv1dx2 * stateDiff[dx2dv10];
    diffs[dv2dv10] = _dv2dx1 * stateDiff[dx1dv10] + _dv2dx2 * stateDiff[dx2dv10];

    diffs[dv1dx10] = _dv1dx1 * stateDiff[dx1dx10] + _dv1dx2 * stateDiff[dx2dx10];
    diffs[dv2dx10] = _dv2dx1 * stateDiff[dx1dx10] + _dv2dx2 * stateDiff[dx2dx10];

    diffs[dv1dv20] = _dv1dx1 * stateDiff[dx1dv20] + _dv1dx2 * stateDiff[dx2dv20];
    diffs[dv2dv20] = _dv2dx1 * stateDiff[dx1dv20] + _dv2dx2 * stateDiff[dx2dv20];

    diffs[dv1dx20] = _dv1dx1 * stateDiff[dx1dx20] + _dv1dx2 * stateDiff[dx2dx20];
    diffs[dv2dx20] = _dv2dx1 * stateDiff[dx1dx20] + _dv2dx2 * stateDiff[dx2dx20];

    // diffs[dv1dk1] = _dv1dk1 + _dv1dx1 * stateDiff[dx1dk1] + _dv1dx2 * stateDiff[dx2dk1];
    // diffs[dv2dk1] = _dv2dx1 * stateDiff[dx1dk1] + _dv2dx2 * stateDiff[dx2dk1];

    // diffs[dv1dk2] = _dv1dk2 + _dv1dx1 * stateDiff[dx1dk2] + _dv1dx2 * stateDiff[dx2dk2];
    // diffs[dv2dk2] = _dv2dk2 + _dv2dx1 * stateDiff[dx1dk2] + _dv2dx2 * stateDiff[dx2dk2];

    // diffs[dv1dl1] = _dv1dl1 + _dv1dx1 * stateDiff[dx1dl1] + _dv1dx2 * stateDiff[dx2dl1];
    // diffs[dv2dl1] = _dv2dx1 * stateDiff[dx1dl1] + _dv2dx2 * stateDiff[dx2dl1];

    // diffs[dv1dl2] = _dv1dl2 + _dv1dx1 * stateDiff[dx1dl2] + _dv1dx2 * stateDiff[dx2dl2];
    // diffs[dv2dl2] = _dv2dl2 + _dv2dx1 * stateDiff[dx1dl2] + _dv2dx2 * stateDiff[dx2dl2];

    diffs[dv1dm1] = _dv1dm1 + _dv1dx1 * stateDiff[dx1dm1] + _dv1dx2 * stateDiff[dx2dm1];
    diffs[dv2dm1] = _dv2dx1 * stateDiff[dx1dm1] + _dv2dx2 * stateDiff[dx2dm1];

    diffs[dv1dm2] = _dv1dx1 * stateDiff[dx1dm2] + _dv1dx2 * stateDiff[dx2dm2];
    diffs[dv2dm2] = _dv2dm2 + _dv2dx1 * stateDiff[dx1dm2] + _dv2dx2 * stateDiff[dx2dm2];

}


void solveStatePlusJacODEs( 
    std::vector<double> const & params,
    double odeIntegrationStep,
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
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0
    };
    double t = 0;
    
    // double debugDiffs[STATE_SIZE + STATE_SIZE*PARAM_SIZE];
    for (auto tStep = 0; tStep < sampleSteps; ++tStep) {

        double* p = curSol.data() + STATE_SIZE;
        for(auto i = 0; i < STATE_SIZE; ++i) {
            state.matrix[tStep * STATE_SIZE + i][0] = curSol[i];
            diffStateByParams.matrix[tStep * STATE_SIZE + i].assign(
                    p,  p+PARAM_SIZE
                );
            p += PARAM_SIZE;
        }
        
        dormandPrince(
            t,
            curSol.data(),
            odeIntegrationStep,
            t,
            curSol.data(),
            params.data(),
            calculateDiffs,
            STATE_SIZE + STATE_SIZE*PARAM_SIZE
        );
    }
}

// Matrix residual();

#endif