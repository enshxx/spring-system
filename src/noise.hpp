#include <random>
#include <cmath>
#include "matrix.hpp"


void generate_noise(Matrix &matrix, double sigma = 1.0) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<double> dist(0.0, sigma);

    for (size_t i = 0; i < matrix.rows(); ++i)
        for (size_t j = 0; j < matrix.cols(); ++j)
            matrix(i, j) = dist(gen) + matrix(i, j);
}
