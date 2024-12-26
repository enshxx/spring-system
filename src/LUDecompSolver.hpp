#ifndef LUDECOMPSOLVER_H
#define LUDECOMPSOLVER_H
#include <iostream>
#include <vector>
#include "matrix.hpp"
// Function to perform LU decomposition
void luDecomposition(Matrix &matrix, int n,
                     Matrix &L,
                     Matrix &U) {
  L = Matrix(n, n);
  U = Matrix(n, n);

  for (int i = 0; i < n; i++) {
    L.matrix[i][i] = 1.0; // diagonal elements of L are 1
    for (int j = i; j < n; j++) {
      U.matrix[i][j] = matrix.matrix[i][j]; // initialize U with original matrix
      for (int k = 0; k < i; k++) {
        U.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
      }
    }
    for (int j = i + 1; j < n; j++) {
      L.matrix[j][i] = matrix.matrix[j][i] / U.matrix[i][i]; // calculate L elements
      for (int k = 0; k < i; k++) {
        L.matrix[j][i] -= L.matrix[j][k] * U.matrix[k][i];
      }
    }
  }
}

// Function to solve linear system using LU decomposition
void solveLinearLU(const Matrix &L,
                       const Matrix &U,
                       const Matrix &b,
                       Matrix &x,
                       int n) {
  Matrix y(n, 1);

  // Forward substitution to solve Ly = b
  for (int i = 0; i < n; i++) {
    y.matrix[i][0] = b.matrix[i][0];
    for (int j = 0; j < i; j++) {
      y.matrix[i][0] -= L.matrix[i][j] * y.matrix[j][0];
    }
    y.matrix[i][0] /= L.matrix[i][i];
  }

  // Backward substitution to solve Ux = y
  x = Matrix(n, 1);
  for (int i = n - 1; i >= 0; i--) {
    x.matrix[i][0] = y.matrix[i][0];
    for (int j = i + 1; j < n; j++) {
      x.matrix[i][0] -= U.matrix[i][j] * x.matrix[j][0];
    }
    x.matrix[i][0] /= U.matrix[i][i];
  }
}

#endif