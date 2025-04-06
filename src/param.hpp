#ifndef PARAM_HPP
#define PARAM_HPP

#include "matrix.hpp"

Matrix create_param(std::initializer_list<double> args) {
    Matrix param(args.size(), 1);

    auto it = args.begin();
    for (size_t i = 0; i < args.size(); ++i) {
        param(i, 0) = *it;
        ++it;
    }

    return param;
}

Matrix create_const_param_for_integrate(double k1, double k2, double m1, double m2, double l1, double l2) {
    return create_param({k1, k2, m1, m2, l1, l2});
}

#endif