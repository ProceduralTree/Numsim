#ifndef PRESSURESOLVERS_H_
#define PRESSURESOLVERS_H_

#include "grid/grid.h"
#include <pde/system.h>
#include <utils/index.h>

struct CGSolver
{
  Grid2D residual;
  Grid2D search_direction;
};

void gauss_seidel_step(PDESystem& system, Index I);

void cg_iteration(PDESystem& system, CGSolver& cg);

#endif // PRESSURESOLVERS_H_
