#include "grid.h"
#include "indexing.h"
#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>

template <Indexing I>
Grid2D<I>::Grid2D(uint16_t x, uint16_t y) : _data(), size_x(x), size_y(y) {

    uint32_t size = std::bit_width(x) + std::bit_width(y) - 2;
    this->_data.resize(1 << size, 0.);
    std::vector<double> &data = this->_data;
    std::fill(data.begin(), data.end(), 0);
};

template <>
double &Grid2D<Indexing::Cartesian>::operator[](uint16_t x, uint16_t y) {
    uint32_t index = x + this->size_x * y;
    return this->_data[index];
}

template <>
const double &Grid2D<Indexing::Cartesian>::operator[](uint16_t x,
                                                      uint16_t y) const {
    uint32_t index = x + this->size_x * y;
    return this->_data[index];
}

template <>
double &Grid2D<Indexing::ZOrder>::operator[](uint16_t x, uint16_t y) {
    uint32_t index = z_order(x, y);
    return this->_data[index];
};

template <>
const double &Grid2D<Indexing::ZOrder>::operator[](uint16_t x,
                                                   uint16_t y) const {
    uint32_t index = z_order(x, y);
    return this->_data[index];
};

template <Indexing I> double &Grid2D<I>::operator[](uint32_t index) {
    return this->_data[index];
};

template <Indexing I>
const double &Grid2D<I>::operator[](uint32_t index) const {
    return this->_data[index];
};

template <Indexing I>
std::ostream &operator<<(std::ostream &os, const Grid2D<I> &obj) {
    os << std::endl;
    os << obj.size_x << "x" << obj.size_y << " Grid2D" << std::endl;
    for (uint16_t i = 0; i < obj.size_x; i++) {
        if (i < 10 || i > obj.size_x - 10) {
            for (uint16_t j = 0; j < obj.size_y; j++) {
                if (j < 10 || j > obj.size_x - 10) {
                    os << std::round(obj[i, j] * 1e5) / 1e5 << "\t";
                } else if (j == 10) {
                    os << "…\t";
                }
            }
            os << std::endl;
        } else if (i == 10) {
            for (uint16_t j = 0; j < obj.size_x; j++) {
                if (j < 10 || j > obj.size_x - 10) {
                    os << "⋮\t";
                }
            }
            os << std::endl;
        }
    }
    return os;
}

bool boundary(uint32_t zindex, uint16_t sx, uint16_t sy) {
    auto [x, y] = decode_z_order(zindex);
    return ((x < sx) && (x > 0)) || ((y < sy) && (y > 0));
}
template <Indexing I>
void print_info(const Grid2D<I> &in, Grid2D<I> &out, uint16_t &z) {
    auto [x, y] = detangle(z);
    auto Idx = indices(z);
    std::cout << "x=" << x << "\t y=" << y << std::endl;
    std::cout << "Top" << Idx.top << "\t" << std::endl;
    std::cout << "Bottom: " << Idx.bottom << "\t" << std::endl;
    std::cout << "Left: " << Idx.left << "\t" << std::endl;
    std::cout << "Right: " << Idx.right << "\t" << std::endl;
    std::cout << "Index: " << z << "\t" << std::endl;
    std::cout << "Elements: " << in.elements() << "\t" << std::endl;
    std::cout << "Out Elements: " << out.elements() << "\t" << std::endl;
    std::cout << "Boundary: " << boundary(z, in.size_x, in.size_y) << "\t"
              << std::endl;
}

template <>
void laplace(const Grid2D<Indexing::ZOrder> &in,
             Grid2D<Indexing::ZOrder> &out) {
    uint32_t N = in.elements();
#pragma omp parallel for
    for (uint32_t z = 1; z < N; z++) {
        if (boundary(z, in.size_x, in.size_y))
            continue;
        auto Idx = indices(z);
        // if (Idx.top > N || Idx.bottom > N || Idx.left > N || Idx.right > N)
        // continue;
        out[z] = (in[Idx.top] + in[Idx.bottom] + in[Idx.left] + in[Idx.right]) -
                 4 * in[z];
    }
}
template <>
void laplace(const Grid2D<Indexing::Cartesian> &in,
             Grid2D<Indexing::Cartesian> &out) {
#pragma omp parallel for
    for (uint16_t j = 0; j < in.size_y; j++) {
        for (uint16_t i = 0; i < in.size_x; i++) {
            if (i == 0 || i == in.size_x - 1 || j == 0 || j == in.size_y - 1)
                continue;
            out[i, j] =
                (in[i + 1, j] + in[i, j + 1] + in[i - 1, j] + in[i, j - 1]) -
                4 * in[i, j];
        }
    }
}
