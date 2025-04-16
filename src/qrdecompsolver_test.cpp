#include "qrdecompsolver.hpp"
#include <iostream>
#include <vector>
#include "matrix.hpp"

int main() {
  Matrix A;
  A.matrix = {
    {4, 12, -16}, {12, 37, -43}, {-16, -43, 98}
  };
  Matrix b(3,1);
  b.matrix = {{10.0}, {20.0}, {30.0}};
  Matrix Q(3, 3), R(3, 3), sol(3,1);
  qrDecompositionSVD(A, Q, R);
  std::cout << "Q matrix" << std::endl;
  printMatrix(Q);
  std::cout << "R matrix " << std::endl;
  printMatrix(R);
  solveLinearQR(Q, R, b, sol, 3);
  std::cout << "solution" << std::endl;
  printMatrix(sol);
  Matrix check(A.rows(), sol.cols());
  multiplyMatrices(A, sol, check);
  std::cout << "expected right hand side" << std::endl;
  printMatrix(b);
  std::cout << "actual right hand side" << std::endl;
  printMatrix(check);
  double tol = 1e-12;
  if (!compareMatrices(check, b, tol)) {
    return 1;
  }
  return 0;
}
