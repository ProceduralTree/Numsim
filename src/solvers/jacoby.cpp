#include "jacoby.h"
#include "grid.h"
#include <memory>
#include <utility>

void jacoby(const Grid2D &b, Grid2D &x, Grid2D &tmp) {

#pragma omp parallel for simd
    for (uint32_t z = 0; z < b.elements(); z++) {
        if (boundary(z, b.elements()))
            continue;

        tmp[z] = 1 / (-4. * x[z]) b[top(z)] -
                 (x[bottom(z)] + x[left(z)] + x[right(z)]);
        std::swap(tmp, x);
    }
}
void jacoby(const Grid2D &b, Grid2D &x) {
    auto tmp = std::make_unique<Grid2D>(b.size_x, b.size_y);
    jacoby(b, x, tmp)
}

void black_red_gauss_seidel(const Grid2D &b, Grid2D &x) {}

void gauss_seidel(PDESystem &system) {
    auto &p = system.p;
    auto &hx = system.hx;
    auto &hy = system.p;
    for (uint16_t j = 1; j < system.size_y - 1; j++) {
        for (uint16_t i = 1; i < system.size_x - 1; i++) {
            p[i, j] = system.b[i, j];
            p[i, j] -= hx * hx * (p[i + 1, j] + p[i - 1, j]);
            p[i, j] -= hy * hy * (p[i, j + 1, j] + p[i, j - 1]);
            p[i, j] /= -2 * hy * hy * p[i, j] - 2 * hx * hx * p[i, j];
        }
    }
}
