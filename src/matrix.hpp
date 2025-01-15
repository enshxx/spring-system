#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>


struct Matrix {
  std::vector<std::vector<double> > matrix;

  static Matrix diag(std::vector<double> &v) {
    Matrix D(v.size(), v.size());
    for (size_t i = 0; i < v.size(); ++i) {
      D.matrix[i][i] = v[i];
    }
    return D;
  }
  static Matrix identity(size_t n) {
    Matrix I(n, n);
    for (size_t i = 0; i < n; ++i) {
      I.matrix[i][i] = 1.0;
    }
    return I;
  }

  Matrix() {}

  Matrix(size_t rows, size_t cols) : 
    matrix(rows, std::vector<double>(cols)) {}
  size_t rows() const {
    return matrix.size();
  }

  void setToZero() {
    for (auto &row : matrix) {
            for (auto &el : row) {
                el = 0.0;
            }
        }
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

// Function to subtract matrix B from matrix A
void subtractMatrices(const Matrix &A,
                      const Matrix &B,
                      Matrix &C) {
  int rowsA = A.rows();
  int colsA = A.cols();

  // Check if matrices can be subtracted
  if (rowsA != B.rows() || colsA != B.cols()) {
    std::cerr << "Error: Matrices cannot be subtracted" << std::endl;
    return;
  }

  // Perform subtraction
  for (size_t i = 0; i < rowsA; i++) {
    for (size_t j = 0; j < colsA; j++) {
      C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
    }
  }
}
// Function to add two matrices
void addMatrices(const Matrix &A,
                const Matrix &B,
                Matrix &C) {
  int rowsA = A.rows();
  int colsA = A.cols();

  // Check if matrices can be added
  if (rowsA != B.rows() || colsA != B.cols()) {
    std::cerr << "Error: Matrices cannot be added" << std::endl;
    return;
  }
  if (rowsA != C.rows() || colsA != C.cols()) {
    std::cerr << "result matrix has incorrect dimensions" << std::endl;
  }

  // Add matrices
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsA; j++) {
      C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
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
  
  if (&C == &A) {
    std::cerr << "Can't use A as result matrix" << std::endl;
    return;
  }
  if (&C == &B) {
    std::cerr << "Can't use B as result matrix" << std::endl;
    return;
  }

  // Check if matrices can be multiplied
  if (colsA != B.rows()) {
    std::cerr << "Error: Matrices cannot be multiplied" << std::endl;
    return;
  }

  if (rowsA != C.rows() || colsB != C.cols()) {
    std::cerr << "result matrix has incorrect dimensions" << std::endl;
    return;
  }

  // Perform matrix multiplication
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsB; j++) {
      C.matrix[i][j] = 0;
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
void transpose(const Matrix &A, Matrix  &B) {
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      B.matrix[j][i] = A.matrix[i][j];
    }
  }
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

// Function to print matrix to file
void printMatrixToFile(const Matrix &mat, const std::string& filename) {
  std::ofstream of(filename, std::ios::out | std::ios::trunc);
  if (!of.good()) {
    of.open(filename, std::ios::out | std::ios::trunc | std::ios::app);
  }
  if (!of) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      of << std::scientific << std::setprecision(6) << std::setw(15) 
      << mat.matrix[i][j] << " ";
    }
    of << std::endl;
  }
  of.close();
}

#endif