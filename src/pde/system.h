#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "grid.h"
#include "indexing.h"
#include <cstdint>

struct Gridsize {
    const double x;
    const double y;
    const double x_squared;
    const double y_squared;

    Gridsize(double x, double y)
        : x(x), y(y), x_squared(x * x), y_squared(y * y) {}
};

struct PDESystem {
    const double Re;
    double dt;
    Grid2D p;
    Grid2D u;
    Grid2D v;
    Grid2D F;
    Grid2D G;
    Grid2D b;
    const uint16_t size_x;
    const uint16_t size_y;
    const Gridsize h;

    PDESystem(double Re, double dt, uint16_t size_x, uint16_t size_y, double hx,
              double hy)
        : Re(Re), dt(dt) , size_x(size_x) , size_y(size_y), p(Grid2D(size_x, size_y)), u(Grid2D(size_x, size_y)),
          v(Grid2D(size_x, size_y)), F(Grid2D(size_x, size_y)),
          G(Grid2D(size_x, size_y)), b(Grid2D(size_x, size_y)),
          h(Gridsize(hx, hy)) {}
};

void print_pde_system(const PDESystem &sys);

#endif // SYSTEM_H_
