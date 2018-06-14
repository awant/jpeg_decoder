#pragma once

#include <vector>
#include <ostream>
#include <cmath>
#include <fftw3.h>


namespace {
    struct Point {
        int x = 0, y = 0;
    };

    struct Direction {
        int x = 1, y = 0;
    };

    void ChangeDirection(const Point &cur_place, Direction *dir, size_t size) {
        if (dir->x == 1 and dir->y == 0) {
            if (cur_place.y > 0) {
                dir->y = -1;
            } else {
                dir->x = -1;
                dir->y = 1;
            }
            return;
        }
        if (dir->x == 0 and dir->y == 1) {
            if (cur_place.x > 0) {
                dir->x = -1;
            } else {
                dir->x = 1;
                dir->y = -1;
            }
            return;
        }
        if (dir->x == 1 and dir->y == -1) {
            if (cur_place.y == 0) {
                dir->y = 0;
            } else if (cur_place.x == (size - 1)) {
                dir->x = 0;
                dir->y = 1;
            }
            return;
        }
        if (dir->x == -1 and dir->y == 1) {
            if (cur_place.y == (size - 1)) {
                dir->x = 1;
                dir->y = 0;
            } else if (cur_place.x == 0) {
                dir->x = 0;
            }
        }
    };
}


template <class T>
class Matrix {
public:
    Matrix(size_t height, size_t width, const T& default_value)
            : height_(height), width_(width),
              buffer_(height, std::vector<T>(width, default_value)) {}

    Matrix(size_t size, const T& default_value)
            : Matrix(size, size, default_value) {}

    Matrix(std::vector<T> array, size_t height, size_t width)
            : Matrix(height, width, T()) {
        if (array.size() != height * width) {
            throw std::runtime_error("Array size doesn't correspond to height, width params");
        }
        size_t idx = 0;
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                buffer_[i][j] = array[idx++];
            }
        }
    }

    T& at(size_t height_coord, size_t width_coord) {
        return buffer_[height_coord][width_coord];
    }

    const T& at(size_t height_coord, size_t width_coord) const {
        return buffer_[height_coord][width_coord];
    }

    size_t GetWidth() const {
        return width_;
    }

    size_t GetHeight() const {
        return height_;
    }

    template <class T2>
    Matrix<T>& Multiply(const Matrix<T2>& rhs) {
        if ((GetWidth() != rhs.GetWidth()) || (GetHeight() != rhs.GetHeight())) {
            throw std::runtime_error("Wrong sizes");
        }
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                buffer_[i][j] *= rhs.at(i, j);
            }
        }
        return *this;
    }

    template <class T2>
    Matrix<T>& Add(const Matrix<T2>& rhs) {
        if ((GetWidth() != rhs.GetWidth()) || (GetHeight() != rhs.GetHeight())) {
            throw std::runtime_error("Wrong sizes");
        }
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                buffer_[i][j] += rhs.at(i, j);
            }
        }
        return *this;
    }

    Matrix<T>& Dot(double coeff) {
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                buffer_[i][j] *= coeff;
            }
        }
        return *this;
    }

    void Fill(T* array) const {
        size_t idx = 0;
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                array[idx++] = buffer_[i][j];
            }
        }
    }

    bool operator==(const Matrix<T>& other) const {
        if ((width_ != other.width_) || (height_ != other.height_)) {
            return false;
        }
        for (int i = 0; i < buffer_.size(); ++i) {
            if (buffer_[i] != other.buffer_[i]) {
                return false;
            }
        }
        return true;
    }

    void Dump(std::ostream& os = std::cout) {
        os << "size: " << height_ << "x" << width_ << "\n";
        os << std::hex;
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                os << buffer_[i][j] << "\t";
            }
            os << "\n";
        }
    }

protected:
    size_t height_, width_;
    std::vector<std::vector<T>> buffer_;
};


template <class T>
class SquareMatrix: public Matrix<T> {
public:
    static SquareMatrix<T> CreateFromZigZag(const std::vector<int>& values, size_t size,
                                            T default_value) {
        SquareMatrix<T> result(size, default_value);
        result.FillZigZag(values);
        return result;
    }

    SquareMatrix(size_t size, const T& default_value)
            : Matrix<T>(size, default_value) {}

    SquareMatrix(std::vector<T> array, size_t size)
            : Matrix<T>(array, size, size) {}

    size_t GetSize() const {
        return this->width_;
    }

    Matrix<T>& FillZigZag(const std::vector<T>& values) {
        if (values.size() > GetSize() * GetSize()) {
            throw std::runtime_error("Matrix can't accommodate result");
        }
        Point pos{0, 0};
        Direction dir{1, -1};
        for (const auto& val : values) {
            this->buffer_[pos.y][pos.x] = val;
            ChangeDirection(pos, &dir, this->height_);
            pos.x += dir.x;
            pos.y += dir.y;
        }
        return *this;
    }
};

//template <class T>
//Matrix<double> IDCT_JPG(Matrix<T>& matrix) {
//    if (matrix.GetWidth() != matrix.GetHeight()) {
//        throw std::runtime_error("Support only square matrices");
//    }
//    size_t size = matrix.GetHeight();
//    Matrix<double> coeffs(size, 1.0);
//    coeffs.at(0, 0) = 0.5;
//    for (size_t i = 0; i < size; ++i) {
//        coeffs.at(0, i) = 1/sqrt(2);
//        coeffs.at(i, 0) = 1/sqrt(2);
//    }
//    coeffs.Multiply(matrix);
//
//    auto* in = new double[size*size];
//    coeffs.Fill(in);
//
//    fftw_plan plan = fftw_plan_r2r_2d(size, size, in, in, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
//    fftw_execute(plan);
//    fftw_destroy_plan(plan);
//
//    Matrix<double> result(in, size, size);
//    delete[] in;
//
//    return result.Dot(0.25);
//}
