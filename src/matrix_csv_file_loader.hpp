#ifndef MATRIX_CSV_FILE_LOADER_HPP
#define MATRIX_CSV_FILE_LOADER_HPP


#include "matrix.hpp"
#include <fstream>
#include <stdexcept>
#include <sstream>


#ifdef PROJECT_SOURCE_DIR
    std::string project_dir = PROJECT_SOURCE_DIR;
#else
    #warning "Макрос PROJECT_SOURCE_DIR не определен"
#endif


class MatrixCSVWriter {

    std::ofstream file;
    std::string filename;

public:
    MatrixCSVWriter(const std::string &filename, std::ios::openmode mode = std::ios::out, bool in_project_dir = true) {
        std::string project_dir;

        if (in_project_dir) {
            #ifdef PROJECT_SOURCE_DIR
                project_dir = PROJECT_SOURCE_DIR;
            #else
                #warning "Макрос PROJECT_SOURCE_DIR не определен"
                project_dir = "";
            #endif
        }
        else
            project_dir = "";

        this->filename = project_dir + "/" + filename;
        file = std::ofstream(this->filename, mode);
    }

    ~MatrixCSVWriter() {
        if (file.is_open())
            file.close();
    }

    MatrixCSVWriter &set_headers(const std::vector<std::string> &headers) {
        if (!file.is_open())
            throw std::runtime_error("Не удалось открыть файл для записи: " + filename);

        for (auto &header : headers)
            file << header << ",";
        file << '\n';

        return *this;
    }

    MatrixCSVWriter &operator<<(const Matrix &matrix) {

        if (!file.is_open())
            throw std::runtime_error("Не удалось открыть файл для записи: " + filename);

        for (size_t i = 0; i < matrix.rows(); ++i) {
            for (size_t j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j < matrix.cols() - 1)
                    file << ",";
            }
            file << '\n';
        }

        return *this;
    }
};


class MatrixCSVRead {

    std::ifstream file;
    std::string filename;

public:

    MatrixCSVRead(const std::string &filename, bool in_project_dir = true) {
        std::string project_dir;

        if (in_project_dir) {
            #ifdef PROJECT_SOURCE_DIR
                project_dir = PROJECT_SOURCE_DIR;
            #else
                #warning "Макрос PROJECT_SOURCE_DIR не определен"
                project_dir = "";
            #endif
        }
        else
            project_dir = "";

        this->filename = project_dir + "/" + filename;
        file = std::ifstream(this->filename);
    }

    ~MatrixCSVRead() {
        if (file.is_open())
            file.close();
    }

    MatrixCSVRead &operator>>(Matrix &matrix) {
        if (!file.is_open())
            throw std::runtime_error("Не удалось открыть файл для чтения: " + filename);

        file.seekg(0, std::ios::beg);

        std::vector<std::vector<long double>> data;
        std::string headers;
        std::getline(file, headers);

        std::string line;

        while (std::getline(file, line)) {
            std::vector<long double> row;
            std::stringstream ss(line);
            std::string value;

            while (std::getline(ss, value, ',')) {
                row.push_back(std::stold(value));
            }

            data.push_back(row);
        }

        if (data.empty())
            throw std::runtime_error("Файл пуст или имеет неправильный формат: " + filename);


        size_t rows = data.size();
        size_t cols = data[0].size();

        matrix = Matrix(rows, cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                matrix(i, j) = data[i][j];
            }
        }

        return *this;
    }
};

#endif