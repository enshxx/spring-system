#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>


struct Matrix {
  std::vector<std::vector<double> > matrix;

  Matrix() {}

  Matrix(size_t rows, size_t cols) : 
    matrix(rows, std::vector<double>(cols)) {}
  size_t rows() const {
    return matrix.size();
  }

  size_t cols() const {
    return matrix.empty() ? 0 : matrix[0].size();
  }
};
// Function to multiply matrix by number
void multiplyMatrixByNumber(Matrix &A, double n) {
 
  for (size_t i = 0; i < A.rows(); ++i) {
    for (size_t j = 0; j < A.cols(); ++j) {
       A.matrix[i][j] *= n;
    }
  }
}
// Function to multiply two matrices
void multiplyMatrices(const Matrix &A,
                      const Matrix &B,
                      Matrix &C) {
  int rowsA = A.rows();
  int colsA = A.cols();
  int colsB = B.cols();

  // Check if matrices can be multiplied
  if (colsA != B.rows()) {
    std::cerr << "Error: Matrices cannot be multiplied" << std::endl;
    return;
  }

  // Initialize result matrix
  C = Matrix(rowsA, colsB);

  // Perform matrix multiplication
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsB; j++) {
      for (int k = 0; k < colsA; k++) {
        C.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
      }
    }
  }
}


// Function to compare two matrices for equality
bool compareMatrices(const Matrix &A, const Matrix &B, double tol) {
  if (A.rows() != B.rows() || A.cols() != B.cols()) {
    std::cout << std::scientific << std::setprecision(15) << "Matrices have different dimensions" << std::endl;
    return false; // Matrices have different dimensions
  }

  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      if (std::abs(A.matrix[i][j] - B.matrix[i][j]) > tol) {
        std::cout << std::scientific << std::setprecision(15) << "Matrices are not equal at row " <<
         i << ", col " << j << " " 
         << A.matrix[i][j] << " !=  "
         <<  B.matrix[i][j] << " diff= "
         << std::abs(A.matrix[i][j] - B.matrix[i][j]) << std::endl;
        return false; // Found an element that is not equal
      }
    }
  }

  return true; // All elements are equal
}
Matrix transpose(const Matrix &A) {
  Matrix B = Matrix(A.cols(), A.rows());
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      B.matrix[j][i] = A.matrix[i][j];
    }
  }
  return B;
}
// Function to print matrix
void printMatrix(const Matrix &mat) {
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      std::cout << std::scientific << mat.matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

#endif