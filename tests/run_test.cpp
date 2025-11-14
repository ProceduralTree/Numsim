
#include <grid/grid.h>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <utils/broadcast.h>

#define ASSERT(condition, message)                               \
  do                                                             \
  {                                                              \
    if (!(condition))                                            \
    {                                                            \
      std::cerr << "Assertion failed: " << message << std::endl; \
      assert(condition);                                         \
    }                                                            \
  } while (0)

void set_one(Index I, Offset O, Grid2D& array)
{

  array[I] += O.x;
  array[I] += O.y;
};

void test_boundary(PDESystem& system)
{
  broadcast_boundary(set_one, system.p.boundary, system.p);
  ASSERT(system.p[system.p.end + Ix] == 1, "system boundary was " << system.p[system.p.end + Ix])

  std::cout << system.p;
}

void test_index()
{
  Index I = { 1, 1 };
  Index Ipx = I + Ix;
  ASSERT(Ipx.x == 2 && Ipx.y == 1, "Plus Failed I.x=" << Ipx.x << " I.y=" << Ipx.y);
  Index Imx = I - Ix;
  ASSERT(Imx.x == 0 && Imx.y == 1, "Minus Failed I.x=" << Imx.x << " I.y=" << Imx.y);
  ASSERT((-5 * Ix).x == -5 && (-5 * Ix).y == 0, "Invert Failed I.x=" << (-5 * Ix).x << " I.y=" << (-5 * Ix).y);
  ASSERT((-Ix).x == -1 && (-Ix).y == 0, "Invert Failed I.x=" << (-Ix).x << " I.y=" << (-Ix).y);
}

int main()
{
  test_index();
};
