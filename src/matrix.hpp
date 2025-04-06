#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iomanip>
#include <fstream>
#include <exception>

class Matrix {

    size_t row, col;
    std::vector<long double> matrix;

public:

    static Matrix diag(const std::vector<double> &v) {
        Matrix D(v.size(), v.size());
        for (size_t i = 0; i < v.size(); ++i)
            D(i, i) = v[i];
        
        return D;
    }

    static Matrix identity(size_t n) {
        Matrix I(n, n);
        for (size_t i = 0; i < n; ++i)
            I(i, i) = 1.0;
        
        return I;
    }

    Matrix(size_t rows, size_t cols) : row(rows), col(cols), matrix(rows * cols, 0) {}

    const long double &operator()(size_t i, size_t j) const {
        if (i >= row || j >= col)
            throw std::out_of_range("Index out of range");

        return matrix[i * col + j];
    }

    long double &operator()(size_t i, size_t j) {
        if (i >= row || j >= col)
            throw std::out_of_range("Index out of range");

        return matrix[i * col + j];
    }

    size_t rows() const { return row; }

    size_t cols() const { return col; }

    void set_zero() {
        for (auto &el : matrix)
            el = 0;
    }

    Matrix operator+(const Matrix &other) const {
        Matrix M(row, col);

        if (row != other.row || col != other.col)
            throw std::invalid_argument("Размер матриц не подходит для сложения");

        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                M(i, j) = (*this)(i, j) + other(i, j);

        return M;
    }

    Matrix &add(const Matrix &other) {
        if (row != other.row || col != other.col)
            throw std::invalid_argument("Размер матриц не подходит для сложения");

        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                (*this)(i, j) += other(i, j);

        return (*this);
    }

    Matrix operator-(const Matrix &other) const {
        Matrix M(row, col);

        if (row != other.row || col != other.col)
            throw std::invalid_argument("Размер матриц не подходит для сложения");

        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                M(i, j) = (*this)(i, j) - other(i, j);

        return M;
    }

    Matrix &sub(const Matrix &other) {
        if (row != other.row || col != other.col)
            throw std::invalid_argument("Размер матриц не подходит для сложения");

        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                (*this)(i, j) -= other(i, j);

        return (*this);
    }

    Matrix operator*(long double n) const {
        Matrix M(row, col);
        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                M(i, j) = (*this)(i, j) * n;
        return M;
    }

    Matrix &mul(long double n) {
        for (auto &el : matrix)
            el *= n;
        return *this;
    }

    Matrix operator*(const Matrix &other) const {
        if (col != other.row)
            throw std::invalid_argument("Размер матриц не подходит для умножения");
        
        Matrix M(row, other.col);

        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < other.col; j++)
                for (size_t k = 0; k < col; k++)
                    M(i, j) += (*this)(i, k) * other(k, j);
        
        return M;
    }

    Matrix &mul(const Matrix &other) {
        if (col != other.row)
            throw std::invalid_argument("Размер матриц не подходит для умножения");

        (*this) = (*this) * other;
        
        return (*this);
    }

    Matrix transpose() const {
        Matrix M(col, row);

        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                M(j, i) = (*this)(i, j);
        
        return M;
    }

    bool equal(const Matrix &other, double sub) const {
        if (row != other.row || col != other.col)
            return false;
        
        for (size_t i = 0; i < row; ++i)
            for (size_t j = 0; j < col; ++j)
                if (abs((*this)(i, j) - other(i, j)) > sub)
                    return false;
        
        return true;
    }

    Matrix &subscribe_row(size_t row, const Matrix &vector) {
        if (this->col != vector.row || vector.col != 1)
            throw std::invalid_argument("Размер вектора не подходит для операции");
        
        if (row >= this->row)
            throw std::out_of_range("Index out of range");
        
        for (size_t i = 0; i < this->col; ++i)
            (*this)(row, i) = vector(i, 0);
        
        return (*this);
    }
};


std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    for (size_t i = 0; i < matrix.rows(); ++i) {
        for (size_t j = 0; j < matrix.cols(); ++j)
            os << std::fixed << std::setprecision(6) << std::setw(12) << matrix(i, j) << " ";
        os << '\n';
    }
    return os;
}

#endif