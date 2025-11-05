#include "output/vtk.h"
#include <grid/grid.h>
#include <grid/indexing.h>
#include <ios>
#include <pde/system.h>

void set_pressure_boundary(PDESystem& system)
{
  auto& p = system.p;
  for (int i = 0; i < system.size_x; i++)
  {
    p[i, 0] = p[i, 1];
  }
  for (int i = 0; i < system.size_x; i++)
  {
    p[i, system.size_y] = p[i, system.size_y - 1];
  }
  for (int j = 0; j < system.size_y; j++)
  {
    p[0, j] = p[1, j];
  }
  for (int j = 0; j < system.size_y; j++)
  {
    p[system.size_x, j] = p[system.size_x - 1, j];
  }
};

void gauss_seidel_step(PDESystem& system, Index I)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = 2 * (1 / h.y_squared) + 2 * (1 / h.x_squared);
  double residual = std::abs(sum_of_neighbours - a_ij * p[I] - system.rhs[I]);
  system.residual = std::max(residual, system.residual);
  p[I] = (sum_of_neighbours - system.rhs[I]) / a_ij;
};

void gauss_seidel(PDESystem& system)
{

  auto& p = system.p;
  auto& h = system.h;
  Grid2D residual(system.size_x + 2, system.size_y + 2);
  std::cout << std::scientific << std::endl;
  std::cout << std::endl;
  std::cout << "Max Pressure: \t" << system.p.max() << "\n";
  std::cout << "Min Pressure: \t" << system.p.min() << "\n";
  std::cout << "Max Velocity x: \t" << system.u.max() << "\n";
  std::cout << "Min Velocity x: \t" << system.u.min() << "\n";
  std::cout << "Max Velocity y: \t" << system.v.max() << "\n";
  std::cout << "Min Velocity y: \t" << system.v.min() << "\n";
  std::cout << "Max RHS: \t" << system.rhs.max() << "\n";
  std::cout << "Min RHS: \t" << system.rhs.min() << "\n";
  std::cout << std::endl;
  // for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
  for (int iter = 0; iter < 10000; iter++)
  {
    std::cout << "Residual : \t" << residual.max() << "\t\r"
              << std::flush;
    set_pressure_boundary(system);
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
      std::cout << "Residual: \t" << residual.max() << "\n";
      break;
    }
  }
  std::cout << std::endl;
}
