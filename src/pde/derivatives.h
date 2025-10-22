#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include "grid.h"
#include "indexing.h"

template <>
constexpr double dx(const &Grid2D field, uint32_t i, uint32_t j, double hx) {
    return 1 / hx * (field[i + 1, j] - field[i, j]);
};

template <>
constexpr double dx_interpolated(const &Grid2D field1, const &Grid2D field2,
                                 uint16_t i, uint16_t j, double hx) {
    return 1 / hx *
           ((field1[i + 1, j] * field2[i + 1, j] +
             field1[i, j] * field2[i, j]) /
                2 -
            (field[i, j] * field2[i, j] + field1[i - 1, j] * field2[i - 1, j]) /
                2);
};

constexpr double ddx(const &Grid2D field, uint_16_t i, uint_16_t j, hx) {
    return 1 / hxsr * (field[i + 1, j] + field[i - 1, j] - 2 * field[i, j]);

#endif // DERIVATIVES_H_
