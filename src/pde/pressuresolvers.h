#ifndef PRESSURESOLVERS_H_
#define PRESSURESOLVERS_H_

#include "grid/grid.h"
#include <pde/system.h>
#include <utils/index.h>

struct CGSolver
{
  Grid2D residual;
  Grid2D search_direction;
  CGSolver(PDESystem& system)
    : residual(system.begin, system.end)
    , search_direction(system.begin, system.end) {
    };
};

struct GaussSeidelSolver
{
  double residual = 0;
};

struct SORSolver
{
  double residual = 0;
};
void gauss_seidel_step(PDESystem& system, Index I);

void cg_iteration(PDESystem& system, CGSolver& cg);

void solve(GaussSeidelSolver gs, PDESystem& system);
void solve(SORSolver gs, PDESystem& system);
void solve(CGSolver& gs, PDESystem& system);
#endif // PRESSURESOLVERS_H_
