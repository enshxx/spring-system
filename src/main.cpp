#include <iostream>
#include <limits>
#include "choleskysolver.hpp"
#include "integrate.hpp"
#include "matrix_csv_file_loader.hpp"
#include "noise.hpp"
#include "param.hpp"
#include "inverse_problem.hpp"

Matrix generate_data(
    double k1, 
    double k2, 
    double m1, 
    double m2, 
    double l1, 
    double l2,
    double h,
    double x1,
    double v1,
    double x2,
    double v2,
    double std_dev,
    size_t iter_num
) {
    Matrix state = create_param({x1, v1, x2, v2});
    Matrix param = create_const_param_for_integrate(k1, k2, m1, m2, l1, l2);

    Matrix states(iter_num, 4);

    for (int i = 0; i < iter_num; ++i) {
        states.subscribe_row(i, state);
        state = dormand_prince(std::function<Matrix(const Matrix&, const Matrix&)>(
                                    compute_derivatives
                                ), state, param, h);
    }

    MatrixCSVWriter("plots/csv_data/integrate_data.csv").set_headers({"x1", "v1", "x2", "v2"}) << states;

    generate_noise(states, std_dev);

    MatrixCSVWriter("plots/csv_data/noise_data.csv").set_headers({"x1", "v1", "x2", "v2"}) << states;

    return states;
}

int main() {
    // Параметры для генерации данных для обратной задачи
    // эти же данные нужны для решения обратной задачи
    double k1 = 2.0, k2 = 10.0, m2 = 50.0, l1 = 6.0, l2 = 1.0;
    double h = 0.01;

    double m1 = 20.0;
    double x1 = 6, v1 = -20.0, x2 = 10.2, v2 = 20.0;
    double std_dev = 0.1;
    size_t iter_num = 1000;


    Matrix states(0, 0);
    // Генерация наблюдений для обратной задачи
    // Полученные данные кладутся в plots/csv_data/integrate_data.csv
    // данные с шумом кладутся в plots/csv_data/noise_data.csv
    generate_data(k1, k2, m1, m2, l1, l2, h, x1, v1, x2, v2, std_dev, iter_num);

    // Чтение наблюдений из файла
    MatrixCSVRead("plots/csv_data/noise_data.csv") >> states;


    // Начальные параметры для решения обратной задачи
    // начальные оценки берутся из наблюдений
    x1 = states(0, 0), v1 = states(0, 1), x2 = states(0, 2), v2 = states(0, 3);
    m1 = (-k1 *(x1 - l1) + k2 * (x2 - x1 - l2));
    
    // Количество итераций Гаусса Ньютона
    iter_num = 1000;

    Matrix state = create_param({x1, v1, x2, v2, m1});
    Matrix param = create_const_param_for_integrate(k1, k2, -std::numeric_limits<double>::infinity(), m2, l1, l2);
    double W = 1 / (std_dev * std_dev);
    InverseProblem solver;
    solver.set_data(states).set_W(W).set_param(param);

    for (int i = 0; i < iter_num; ++i) {
        state = solver.solve(state, h);
        std::cout << state << std::endl;
    }

    return 0;
}
