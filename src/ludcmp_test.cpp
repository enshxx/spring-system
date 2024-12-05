#include "LUDecompSolver.hpp"
#include "matrix.hpp"

int main() {
  int n = 3; // size of matrix
  Matrix matrix(3, 3);
  matrix.matrix = {
      {4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
  Matrix L(3,3), U(3,3);
  Matrix b(3,1);
  b.matrix = {{10}, {20}, {30}}; // right-hand side vector
  Matrix x(n, 1);        // solution vector

  luDecomposition(matrix, n, L, U);

  std::cout << "Original matrix:" << std::endl;
  printMatrix(matrix);

  std::cout << "L matrix:" << std::endl;
  printMatrix(L);

  std::cout << "U matrix:" << std::endl;
  printMatrix(U);

  solveLinearSystem(L, U, b, x, n);

  std::cout << "Solution vector:" << std::endl;
  printMatrix(x);
  std::cout << std::endl;

  Matrix checkMatrix(n, n);
  multiplyMatrices(matrix, x, checkMatrix);
  std::cout << "Check matrix:" << std::endl;
  printMatrix(checkMatrix);
  std::cout << "b matrix:" << std::endl;
  printMatrix(b);
  double tol = 1e-12;
  std::cout << "Comparison tolerance: " << tol << std::endl;
  std::cout << "Machine epsilon: " << std::numeric_limits<double>::epsilon() << std::endl;
  if (!compareMatrices(checkMatrix, b, tol)) {
    std::cerr << "Error: Solution is incorrect" << std::endl;
    return 1;
  }
  return 0;
}
