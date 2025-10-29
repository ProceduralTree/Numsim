#include <grid/grid.h>
#include <grid/indexing.h>
#include <pde/system.h>

void gauss_seidel(PDESystem& system)
{
  auto& p = system.p;
  auto& h = system.h;
  bool converged = false;
  Grid2D residual(system.size_x + 2, system.size_y + 2);
  while (!converged)
  {
    for (uint16_t i = 1; i < system.size_x + 1; i++)
    {
      for (uint16_t j = 1; j < system.size_y + 1; j++)
      {
        double new_p = ((h.x_squared * h.y_squared) / (2 * (h.x_squared + h.y_squared))) * (((p[i - 1, j] + p[i + 1, j]) / h.x_squared) + ((p[i, j - 1] + p[i, j + 1]) / h.y_squared) - system.rhs[i, j]);
        residual[i, j] = std::abs(p[i, j] - new_p);
        p[i, j] = new_p;
      }
    }
    if (residual.max() < 1e-4)
    {
      converged = true;
    }
  }
}
