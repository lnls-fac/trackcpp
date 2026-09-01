// TRACKCPP - Particle tracking code
// Copyright (C) 2015  LNLS Accelerator Physics Group
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef COEFMATRIX_H
#define COEFMATRIX_H

#include <cstddef>
#include <vector>

#ifndef SWIG
template<typename T>
class MatrixRow {
public:
    explicit MatrixRow(T* ptr)
        : ptr_(ptr)
    {}

    T& operator[](size_t j)
    {
        return ptr_[j];
    }

private:
    T* ptr_;
};
#endif


class CoefMatrix {
public:

#ifndef SWIG
    using Row = MatrixRow<double>;
    using ConstRow = MatrixRow<const double>;
#endif

    CoefMatrix() = default;

    CoefMatrix(size_t rows, size_t cols)
        : rows_(rows),
          cols_(cols),
          data_(rows * cols, 0.0)
    {}

    void resize(size_t rows, size_t cols)
    {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols);
    }

    double& operator()(size_t i, size_t j)
    {
        return data_[i * cols_ + j];
    }

    const double& operator()(size_t i, size_t j) const
    {
        return data_[i * cols_ + j];
    }

#ifndef SWIG
    Row operator[](size_t i)
    {
        return Row(data_.data() + i * cols_);
    }

    ConstRow operator[](size_t i) const
    {
        return ConstRow(data_.data() + i * cols_);
    }
#endif

    size_t rows() const
    {
        return rows_;
    }

    size_t cols() const
    {
        return cols_;
    }

    size_t size() const
    {
        return data_.size();
    }

    bool empty() const
    {
        return data_.empty();
    }

    double* data()
    {
        return data_.data();
    }

    const double* data() const
    {
        return data_.data();
    }

private:

    size_t rows_ = 0;
    size_t cols_ = 0;
    std::vector<double> data_;
};

#endif