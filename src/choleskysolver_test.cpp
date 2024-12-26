#include "matrix.hpp"
#include "linearsystemsolver.hpp"

int main() {
    Matrix A(3, 3);
    A.matrix = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    Matrix b(3, 1);
    b.matrix = {{4}, {-10}, {2}};
    Matrix x(3, 1);
    Matrix L(A.rows(), A.rows());
    solveLinearChol(A, L, b, x);
    Matrix realB(3, 1);
    multiplyMatrices(A, x, realB);
    if (!compareMatrices(b, realB, 1e-11)) {
        return 1;
    }
    return 0;
}