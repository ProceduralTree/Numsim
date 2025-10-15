#include "grid.h"
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <iostream>

void hello() { std::cout << "Hello" << std::endl; }

Grid2D::Grid2D(uint32_t x, uint32_t y) : _data(), size_x(x), size_y(y) {
    uint64_t size = std::bit_width(x) + std::bit_width(y);
    this->_data.reserve(1 << size);
    std::vector<double> &data = this->_data;
    std::fill(data.begin(), data.end(), 0);
};
double &Grid2D::operator[](uint32_t x, uint32_t y) {
    uint64_t index = interleave(x, y);
    return this->_data[index];
};

const double &Grid2D::operator[](uint32_t x, uint32_t y) const {
    uint64_t index = interleave(x, y);
    return this->_data[index];
};
uint64_t interleave(uint32_t x, uint32_t y) {
    uint64_t maskx = 0x5555555555555555;
    uint64_t masky = 0xAAAAAAAAAAAAAAAA;
    uint64_t xbits = _pdep_u64(x, maskx);
    uint64_t ybits = _pdep_u64(y, masky);
    return xbits | ybits;
}

std::ostream &operator<<(std::ostream &os, const Grid2D &obj) {
    os << std::endl;
    os << obj.size_x << "x" << obj.size_y << " Grid2D" << std::endl;
    for (uint32_t i = 0; i < obj.size_x; i++) {
        if (i < 10 || i > obj.size_x - 10) {
            for (uint32_t j = 0; j < obj.size_y; j++) {
                if (j < 10 || j > obj.size_x - 10) {
                    os << std::round(obj[i, j] * 1e5) / 1e5 << "\t";
                } else if (j == 10) {
                    os << "…\t";
                }
            }
            os << std::endl;
        } else if (i == 10) {
            for (uint32_t j = 0; j < obj.size_x; j++) {
                if (j < 10 || j > obj.size_x - 10) {
                    os << "⋮\t";
                }
            }
            os << std::endl;
        }
    }
    return os;
}

void laplace(const Grid2D &in, Grid2D &out) {
    assert(in.size_x == out.size_x);
    assert(in.size_y == out.size_y);
    for (uint32_t i = 1; i < in.size_x - 1; i++) {
        for (uint32_t j = 1; j < in.size_y - 1; j++) {
            out[i, j] = in[i + 1, j] + in[i, j + 1] + in[i - 1, j] +
                        in[i, j - 1] - 4 * in[i, j];
        }
    }
}
