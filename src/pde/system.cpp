#include "system.h"
#include "derivatives.h"
#include "utils/broadcast.h"
#include <cstdint>

void calculate_FG(PDESystem &system, uint16_t i, uint16_t j) {
    auto &u = system.u;
    auto &v = system.v;
    auto &h = system.h;
    system.F[i, j] =
        u[i, j] +
        system.dt *
            (1 / system.Re * (ddx(u, i, j, h) + ddy(u, i, j, h)) -
             dx_interpolated(u, u, i, j, h) - dx_interpolated(u, v, i, j, h));
    system.G[i, j] =
        v[i, j] +
        system.dt *
            (1 / system.Re * (ddy(v, i, j, h) + ddx(v, i, j, h)) -
             dy_interpolated(v, v, i, j, h) - dy_interpolated(u, v, i, j, h));
}

void update_uv(PDESystem &system, uint16_t i, uint16_t j) {
    auto &p = system.p;
    auto &h = system.h;
    system.u[i, j] = system.F[i, j] - system.dt * dx(p, i, j, h);
    system.v[i, j] = system.G[i, j] - system.dt * dy(p, i, j, h);
};

void solve_pressure(PDESystem &system, uint16_t i, uint16_t j) {
    std::cout << "Hello World" << std::endl;
};

void timestep(PDESystem system) {
    broadcast(calculate_FG, system);
    broadcast(solve_pressure, system);
    broadcast(update_uv, system);
    //
}
void print_pde_system(const PDESystem &sys) {
    printf("╔═══════════════════════════════════════════════╗\n");
    printf("║              PDE System Summary               ║\n");
    printf("╚═══════════════════════════════════════════════╝\n");

    printf("Reynolds number (Re): %.6f\n", sys.Re);
    printf("Time step (dt):       %.6f\n", sys.dt);
    printf("Grid size (x, y):     %u x %u\n", sys.size_x, sys.size_y);
    printf("Grid spacing (dx, dy): %.6f, %.6f\n", sys.h.x, sys.h.y);
    printf("\n");
    printf("───────────────────────────────────────────────\n");
};
