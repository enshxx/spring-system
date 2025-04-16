#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include "matrix.hpp"
#include <iostream>
#include <vector>
#include <functional>

Matrix compute_derivatives_into_deviations(const Matrix& state, const Matrix& const_state) {
    Matrix derivatives(24, 1);
    /*
    [
        dx1/dt, dv1/dt, dx2/dt, dv2/dt, dx1/dx1, dv1/dx1, dx2/dx1, dv2/dx1, 
        dx1/dv1, dv1/dv1, dx2/dv1, dv2/dv1, dx1/dx2, dv1/dx2, dx2/dx2, dv2/dx2,
        dx1/dv2, dv1/dv2, dx2/dv2, dv2/dv2, dx1/dm1, dv1/dm1, dx2/dm1, dv2/dm1
    ]
    */

    double k1 = const_state(0, 0);
    double k2 = const_state(1, 0);
    double m1 = const_state(2, 0);
    double m2 = const_state(3, 0);
    double l1 = const_state(4, 0);
    double l2 = const_state(5, 0);

    double dx1dt = state(0, 0);
    double dv1dt = state(1, 0);
    double dx2dt = state(2, 0);
    double dv2dt = state(3, 0);
    double dx1dx1 = state(4, 0);
    double dv1dx1 = state(5, 0);
    double dx2dx1 = state(6, 0);
    double dv2dx1 = state(7, 0);
    double dx1dv1 = state(8, 0);
    double dv1dv1 = state(9, 0);
    double dx2dv1 = state(10, 0);
    double dv2dv1 = state(11, 0);
    double dx1dx2 = state(12, 0);
    double dv1dx2 = state(13, 0);
    double dx2dx2 = state(14, 0);
    double dv2dx2 = state(15, 0);
    double dx1dv2 = state(16, 0);
    double dv1dv2 = state(17, 0);
    double dx2dv2 = state(18, 0);
    double dv2dv2 = state(19, 0);
    double dx1dm1 = state(20, 0);
    double dv1dm1 = state(21, 0);
    double dx2dm1 = state(22, 0);
    double dv2dm1 = state(23, 0);


    double _dv1dx1 = -(k1 / m1 + k2 / m1);
    double _dv1dx2 = k2 / m1;
    double _dv1dm1 = k1 * (dx1dt - l1) / (m1 * m1) - k2 * (dx2dt - dx1dt - l2) / (m1 * m1);

    double _dv2dx1 = k2 / m2;
    double _dv2dx2 = -_dv2dx1;
    

    derivatives(0, 0) = dv1dt;
    derivatives(1, 0) = -k1 / m1 * (dx1dt - l1) + k2 / m1 * (dx2dt - dx1dt - l2);
    derivatives(2, 0) = dv2dt;
    derivatives(3, 0) = -k2 / m2 * (dx2dt - dx1dt - l2);

    derivatives(4, 0) = dv1dx1;
    derivatives(5, 0) = _dv1dx1 + _dv1dx2 * dx2dx1;
    derivatives(6, 0) = dv2dx1;
    derivatives(7, 0) = _dv2dx1 + _dv2dx2 * dx2dx1;

    derivatives(8, 0) = dv1dv1;
    derivatives(9, 0) = _dv1dx1 * dx1dv1 + _dv1dx2 * dx2dv1;
    derivatives(10, 0) = dv2dv1;
    derivatives(11, 0) = _dv2dx1 * dx1dv1 + _dv2dx2 * dx2dv1;

    derivatives(12, 0) = dv1dx2;
    derivatives(13, 0) = _dv1dx2 + _dv1dx1 * dx1dx2;
    derivatives(14, 0) = dv2dx2;
    derivatives(15, 0) = _dv2dx2 + _dv2dx1 * dx1dx2;

    derivatives(16, 0) = dv1dv2;
    derivatives(17, 0) = _dv1dx1 * dx1dv2 + _dv1dx2 * dx2dv2;
    derivatives(18, 0) = dv2dv2;
    derivatives(19, 0) = _dv2dx1 * dx1dv2 + _dv2dx2 * dx2dv2;
    
    derivatives(20, 0) = dv1dm1;
    derivatives(21, 0) = _dv1dm1 + _dv1dx1 * dx1dm1 + _dv1dx2 * dx2dm1;
    derivatives(22, 0) = dv2dm1;
    derivatives(23, 0) = _dv2dx1 * dx1dm1 + _dv2dx2 * dx2dm1;

    return derivatives;
}

Matrix compute_derivatives(const Matrix& state, const Matrix& const_state) {
    Matrix derivatives(4, 1); // [dx1/dt, dv1/dt, dx2/dt, dv2/dt]

    double x1 = state(0, 0);
    double v1 = state(1, 0);
    double x2 = state(2, 0);
    double v2 = state(3, 0);

    double k1 = const_state(0, 0);
    double k2 = const_state(1, 0);
    double m1 = const_state(2, 0);
    double m2 = const_state(3, 0);
    double l1 = const_state(4, 0);
    double l2 = const_state(5, 0);

    derivatives(0, 0) = v1;
    derivatives(1, 0) = -k1 / m1 * (x1 - l1) + k2 / m1 * (x2 - x1 - l2);
    derivatives(2, 0) = v2;
    derivatives(3, 0) = -k2 / m2 * (x2 - x1 - l2);

    return derivatives;
}

Matrix dormand_prince(std::function<Matrix(const Matrix&, const Matrix&)> derivative, const Matrix& state, const Matrix& const_state, double h) {
    const size_t stages = 13;
    const double a[stages][stages] = {
        {0},
        {1.0 / 18.0},
        {1.0 / 48.0, 1.0 / 16.0},
        {1.0 / 32.0, 0, 3.0 / 32.0},
        {5.0 / 16.0, 0, -75.0 / 64.0, 75.0 / 64.0},
        {3.0 / 80.0, 0, 0, 3.0 / 16.0, 3.0 / 20.0},
        {29443841.0 / 614563906.0, 0, 0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0},
        {16016141.0 / 946692911.0, 0, 0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0},
        {39632708.0 / 573591083.0, 0, 0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0},
        {246121993.0 / 1340847787.0, 0, 0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0},
        {-1028468189.0 / 846180014.0, 0, 0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0},
        {185892177.0 / 718116043.0, 0, 0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0},
        {403863854.0 / 491063109.0, 0, 0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0}
    };

    const double b[stages] = {
        35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0
    };

    std::vector<Matrix> k(stages, Matrix(state.rows(), 1));

    for (size_t i = 0; i < stages; ++i) {
        Matrix sum(state.rows(), 1);
        for (size_t j = 0; j < i; ++j)
            sum.add(k[j] * a[i][j]);
        
        k[i] = derivative(state + sum * h, const_state);
    }

    Matrix newState = state;
    for (size_t i = 0; i < stages; ++i)
        newState.add(k[i] * (b[i] * h));
    

    return newState;
}

#endif