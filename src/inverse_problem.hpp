#ifndef INVERSE_PROBLEM_HPP
#define INVERSE_PROBLEM_HPP

#include "integrate.hpp"
#include "choleskysolver.hpp"
#include <vector>
#include <cmath>
#include <functional>


class InverseProblem {

    Matrix data;
    double W;
    Matrix param;

public:

    InverseProblem() : data(Matrix(4, 0)), W(1), param(Matrix(6, 0)) {}

    InverseProblem &set_data(const Matrix &data) {
        this->data = data;
        return *this;
    }
    InverseProblem &set_W(double W) {
        this->W = W;
        return *this;
    }
    InverseProblem &set_param(const Matrix &param) {
        this->param = param;
        return *this;
    }

    Matrix solve(const Matrix& state, double h) {
        Matrix AtWA(5, 5);
        Matrix AtWr(5, 1);
        Matrix tmp(24, 1);
        Matrix partial_derivatives(5, 1);
        /*
        [
            dx1/dt, dv1/dt, dx2/dt, dv2/dt, dx1/dx1, dv1/dx1, dx2/dx1, dv2/dx1, 
            dx1/dv1, dv1/dv1, dx2/dv1, dv2/dv1, dx1/dx2, dv1/dx2, dx2/dx2, dv2/dx2,
            dx1/dv2, dv1/dv2, dx2/dv2, dv2/dv2, dx1/dm1, dv1/dm1, dx2/dm1, dv2/dm1
        ]
        */

        tmp(0, 0) = state(0, 0);
        tmp(1, 0) = state(1, 0);
        tmp(2, 0) = state(2, 0);
        tmp(3, 0) = state(3, 0);
        tmp(4, 0) = 1;
        tmp(9, 0) = 1;
        tmp(14, 0) = 1;
        tmp(19, 0) = 1;

        param(2, 0) = state(4, 0);


        for (size_t i = 0; i < data.rows(); ++i) {
            for (size_t j = 0; j < 4; ++j) {
                for (size_t k = 0; k < 5; ++k)
                    partial_derivatives(k, 0) = tmp(j + 4 * k + 4, 0);

                AtWA.add((partial_derivatives * partial_derivatives.transpose()).mul(W));
                AtWr.add(partial_derivatives.mul(data(i, j) - tmp(j, 0)).mul(W));
            }

            tmp = dormand_prince(
                std::function<Matrix(const Matrix&, const Matrix&)>(
                    compute_derivatives_into_deviations
                ),
                tmp, param, h
            );
        }


        LinearSolver solver(AtWA, AtWA * state + AtWr);
        Matrix new_state = solver.solve();
        return new_state;
    }

};

#endif