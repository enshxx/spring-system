#ifndef LINEARSYSTEMSOLVER_H
#define LINEARSYSTEMSOLVER_H
#include "matrix.hpp"
#include "LUDecompSolver.hpp"
#include "qrdecompsolver.hpp"
#include "choleskysolver.hpp"

void solveLinearSystemLU(Matrix &A, Matrix &b, Matrix &x) {
    Matrix L, U;
    luDecomposition(A, A.rows(), L, U);
    solveLinearLU(L, U, b, x, A.rows());
}

void solveLinearSystemQR(Matrix &A, Matrix &b, Matrix &x) {
    Matrix Q, R;
    qrDecomposition(A, Q, R);
    solveLinearQR(Q, R, b, x, A.rows());
}

void solveLinearChol(Matrix &A, Matrix &L, Matrix &b, Matrix &x) {
    if (choleskyDecomposition(A, L)) {
        solveLinearCholesky(L, b, x, L.cols());
    }
}

#endif