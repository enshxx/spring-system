#ifndef     CHOLESKYSOLVER_H
#define     CHOLESKYSOLVER_H
#include <cmath>
#include "matrix.hpp"


// Function to perform Cholesky decomposition
bool choleskyDecomposition(const Matrix &matrix, Matrix &L) {
    int n = matrix.rows();
    if (matrix.cols() != n) {
        std::cerr << "Error: Matrix must be square for Cholesky decomposition." << std::endl;
        return false;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            for (int k = 0; k < j; k++)
                sum += L.matrix[i][k] * L.matrix[j][k];

            if (i == j)
                L.matrix[i][j] = std::sqrt(matrix.matrix[i][i] - sum);
            else
                L.matrix[i][j] = (matrix.matrix[i][j] - sum) / L.matrix[j][j];
        }

        if (L.matrix[i][i] <= 0.0) {
            std::cerr << "Error: Matrix is not positive definite." << std::endl;
            return false;
        }
    }

    return true;
}

// Function to solve linear system using Cholesky decomposition
void solveLinearCholesky(const Matrix &L,
                         const Matrix &b,
                         Matrix &x,
                         int n) {
    Matrix y(n, 1);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L.matrix[i][j] * y.matrix[j][0];
        }
        y.matrix[i][0] = (b.matrix[i][0] - sum) / L.matrix[i][i];
    }

    x = Matrix(n, 1);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L.matrix[j][i] * x.matrix[j][0];
        }
        x.matrix[i][0] = (y.matrix[i][0] - sum) / L.matrix[i][i];
    }
}

#endif