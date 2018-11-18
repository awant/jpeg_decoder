#pragma once

#include <vector>
#include <ostream>
#include <cmath>
#include <fftw3.h>
#include <iostream>

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
            } else if (cur_place.x == static_cast<int>(size - 1)) {
                dir->x = 0;
                dir->y = 1;
            }
            return;
        }
        if (dir->x == -1 and dir->y == 1) {
            if (cur_place.y == static_cast<int>(size - 1)) {
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
    const double epsilon = 0.001;
public:
    Matrix() { }

    Matrix(size_t height, size_t width, const T& default_value)
            : height_(height), width_(width),
              buffer_(height, std::vector<T>(width, default_value)) {}

    Matrix(size_t size, const T& default_value)
            : Matrix(size, size, default_value) {}

    Matrix(size_t height, size_t width, const std::vector<T>& array)
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

    template <class T2>
    explicit Matrix(const Matrix<T2>& matrix): Matrix(matrix.GetHeight(), matrix.GetWidth(), T()) {
        for (int y = 0; y < height_; ++y) {
            for (int x = 0; x < width_; ++x) {
                buffer_[y][x] = matrix.at(y, x);
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

    void Map(const Point& upper_left_corner, const Point& lower_right_corner, const Matrix<T>& rhs) {
        int width = lower_right_corner.x - upper_left_corner.x;
        int height = lower_right_corner.y - upper_left_corner.y;

        assert(rhs.GetWidth() == width);
        assert(rhs.GetHeight() == height);
        assert(lower_right_corner.y <= static_cast<int>(height_));
        assert(lower_right_corner.x <= static_cast<int>(width_));
        // map values
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                buffer_[upper_left_corner.y+y][upper_left_corner.x+x] = rhs.at(y, x);
            }
        }
    }

    bool operator==(const Matrix<T>& rhs) const {
        if ((width_ != rhs.GetWidth()) || (height_ != rhs.GetHeight())) {
            return false;
        }
        for (int i = 0; i < height_; ++i) {
            for (int j = 0; j < width_; j++) {
                if (abs(buffer_[i][j] - rhs.buffer_[i][j]) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    void Dump(std::ostream& os = std::cout) const {
        os << "size: " << height_ << "x" << width_ << "\n";
        os << std::dec;
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
    static SquareMatrix<T> CreateFromZigZag(size_t size, const std::vector<int>& values,
                                            T default_value) {
        SquareMatrix<T> result(size, default_value);
        result.FillZigZag(values);
        return result;
    }

    SquareMatrix() { }

    SquareMatrix(size_t size, const T& default_value)
            : Matrix<T>(size, default_value) {}

    SquareMatrix(size_t size, const std::vector<T>& array)
            : Matrix<T>(size, size, array) {}

    template<class T2>
    explicit SquareMatrix(const SquareMatrix<T2>& rhs)
            : Matrix<T>(rhs) {}

    size_t GetSize() const {
        return this->width_;
    }

    Matrix<T>& FillZigZag(const std::vector<int>& values) {
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


template <class T>
class MatrixTransformer {
public:
    static const int TransformSize = 8;
    const double coeff1 = 0.5 / TransformSize;
    const double coeff2 = std::sqrt(2);

    MatrixTransformer() {
        plan_ = fftw_plan_r2r_2d(TransformSize, TransformSize, input_buffer_,
                output_buffer_, FFTW_REDFT01, FFTW_REDFT01, 0);
    };

    ~MatrixTransformer() {
        fftw_destroy_plan(plan_);
    }

    void MakeIDCTransform(SquareMatrix<T>* matrix) {
        if (matrix->GetSize() != TransformSize) {
            throw std::runtime_error("Support only 8x8 matrices");
        }

        // fill input buffer
        for (size_t i = 0; i < TransformSize; ++i) {
            for (size_t j = 0; j < TransformSize; ++j) {
                input_buffer_[i*TransformSize+j] = coeff1 * matrix->at(i, j);
                if (!i) {
                    input_buffer_[i*TransformSize+j] *= coeff2;
                }
                if (!j) {
                    input_buffer_[i*TransformSize+j] *= coeff2;
                }
            }
        }

        fftw_execute(plan_);

        // fill matrix from output buffer
        for (size_t i = 0; i < TransformSize; ++i) {
            for (size_t j = 0; j < TransformSize; ++j) {
                matrix->at(i, j) = output_buffer_[i*TransformSize+j];
            }
        }
    }

private:
    double input_buffer_[TransformSize*TransformSize];
    double output_buffer_[TransformSize*TransformSize];
    fftw_plan plan_;
};
