#include "system.h"
#include "broadcast.h"
#include "derivatives.h"

void timestep(PDESystem system) {
    broadcast(calculate_FG, system);
    broadcast(solve_pressure, system);
    broadcast(update_uv, system);
    //
}

void calculate_FG(PDESystem &system, uint16_t i, uint16_t j) {
    auto &u = system.u;
    auto &v = system.v;
    system.F[i, j] =
        u[i, j] +
        system.dt * (1 / system.Re * (ddx(u, i, j) + ddy(u, i, j)) -
                     dx_interpolated(u, u, i, j) - dx_interpolated(u, v, i, j));
    system.G[i, j] =
        v[i, j] +
        system.dt * (1 / system.Re * (ddy(v, i, j) + ddx(v, i, j)) -
                     dy_interpolated(v, v, i, j) - dy_interpolated(u, v, i, j));
}

void update_uv(PDESystem system, uint16_t i, uint16_t j) {
    auto &p = system.p;
    system.u[i, j] = system.F[i, j] - system.dt * dx(p, i, j);
    system.v[i, j] = system.G[i, j] - system.dt * dy(p, i, j);
};
