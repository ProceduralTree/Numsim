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
  // Grid2D residual;
  // GaussSeidelSolver(PDESystem& system)
  //   : residual(system.begin, system.end) { };
};

struct SORSolver
{
  // Grid2D residual;
  // SORSolver(PDESystem& system)
  //   : residual(system.begin, system.end) { };
};
// struct BlackRedSolver
//{
//   Grid2D residual;
//   BlackRedSolver(PDESystem& system)
//     : residual(system.begin, system.end) { };
// };

void cg_iteration(PDESystem& system, CGSolver& cg);

void solve(GaussSeidelSolver& S, PDESystem& system);
void solve(SORSolver& S, PDESystem& system);
void solve(CGSolver& S, PDESystem& system);
#endif // PRESSURESOLVERS_H_
