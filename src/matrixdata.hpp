#ifndef MATRIXDATA_H
#define MATRIXDATA_H
#include "matrix.hpp"

#define STATE_SIZE 4
#define PARAM_SIZE 10
struct MatrixData {
    Matrix FdiffState;
    Matrix FdiffParam;
    Matrix stateDiffThetaDot;
    Matrix tmp;
    Matrix stateDiffThetaMat;
    MatrixData() : 
        FdiffState(STATE_SIZE, STATE_SIZE),
        FdiffParam(STATE_SIZE, PARAM_SIZE),
        stateDiffThetaDot(STATE_SIZE, PARAM_SIZE),
        tmp(STATE_SIZE, PARAM_SIZE),
        stateDiffThetaMat(STATE_SIZE, PARAM_SIZE)
    {}
    
    void setToZero() {
        FdiffState.setToZero();
        FdiffParam.setToZero();
        stateDiffThetaDot.setToZero();
        tmp.setToZero();
        stateDiffThetaMat.setToZero();

    }

};

#endif