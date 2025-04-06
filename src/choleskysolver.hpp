#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include <cmath>
#include "matrix.hpp"
#include <execution>


class LinearSolver {

    Matrix L, b;

    Matrix cholesky_decomposition(const Matrix &matrix) {
        size_t n = matrix.rows();
        if (n != matrix.cols())
            throw std::invalid_argument("Матрица должна быть квадратной");

        Matrix L(n, n);

        for (size_t k = 0; k < n; ++k) {
            double sum = 0.0;
            for (size_t p = 0; p < k; ++p)
                sum += L(k, p) * L(k, p);
            
            L(k, k) = matrix(k, k) - sum;

            // if (L(k, k) <= 0)
            //     throw std::invalid_argument("Матрица не положительно определена");

            L(k, k) = std::sqrt(L(k, k));

            for (size_t i = k + 1; i < n; ++i) {
                sum = 0.0;
                for (size_t p = 0; p < k; ++p)
                    sum += L(i, p) * L(k, p);
                L(i, k) = (matrix(i, k) - sum) / L(k, k);
            }
        }

        return L;
    }

public:

    LinearSolver() : L(Matrix::identity(0)), b(Matrix::identity(0)) {}

    LinearSolver(const Matrix &matrix, const Matrix &vector) : L(Matrix::identity(0)), b(vector) {
        set_matrix(matrix);
    }

    void set_matrix(const Matrix &matrix) {
        try {
            L = cholesky_decomposition(matrix);
        }
        catch (const std::invalid_argument &e) {
            throw std::invalid_argument("Матрица должна быть квадратной");
        }
    }

    void set_vector(const Matrix &vector) { b = vector; }

    Matrix solve() {
        size_t n = L.rows();
        if (b.rows() != n || b.cols() != 1)
            throw std::invalid_argument("Вектор должен быть столбцом с размером, равным размеру матрицы");

        Matrix y(n, 1);
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j)
                sum += L(i, j) * y(j, 0);
            y(i, 0) = (b(i, 0) - sum) / L(i, i);
        }

        Matrix x(n, 1);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (size_t j = i + 1; j < n; ++j)
                sum += L(j, i) * x(j, 0);
            x(i, 0) = (y(i, 0) - sum) / L(i, i);
        }

        return x;
    }

};

#endif