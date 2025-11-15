#include "utils/Logger.h"
#include "utils/broadcast.h"
#include "utils/settings.h"
#include <csignal>
#include <grid/grid.h>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <utils/profiler.h>

#define ASSERT(condition, message)                               \
  do                                                             \
  {                                                              \
    if (!(condition))                                            \
    {                                                            \
      std::cerr << "Assertion failed: " << message << std::endl; \
      assert(condition);                                         \
    }                                                            \
  } while (0)

void handle_sigint(int)
{
  std::cout << "\n[Profiler] Ctrl+C detected, printing results...\n";
  Profiler::instance().print();
  std::_Exit(0); // exit immediately
}

void set_one(Index I, Offset O, Grid2D& array)
{

  array[I] += O.x;
  array[I] += O.y;
};

void test_boundary(PDESystem& system)
{
  broadcast_boundary(set_one, system.p.boundary, system.p);

  // broadcast_x_boundary(
  //   [&](PDESystem& s, Index I, Offset o) {
  //     s.v[I + o] = o.y;
  //   },
  //   system, system.v);
  // std::cout << system.v;
  // std::cout << "\n";
  // broadcast_boundary(
  //   [&](PDESystem& s, Index I, Offset o) {
  //     s.p[I + o] = o.x + o.y;
  //   },
  //   system, system.p);
  std::cout << system.p;
  std::cout << "\n";
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

auto main(int argc, char* argv[]) -> int
{
  std::signal(SIGINT, handle_sigint);
  LOG::Init(LOG::LoggerType::STDOUT);
  if (argc < 2)
  {
    LOG::Warning("missing file name");
    return -1;
  }
  test_index();
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (!Settings::loadFromFile(argv[1]))
  {
    LOG::Warning("couldn't parse settings file");
    return -1;
  }
  Settings::get().printSettings();

  PDESystem test_system = PDESystem(Settings::get());

  print_pde_system(test_system);

  // test_boundary(test_system);
  double time = 0;
  while (time < test_system.settings.endTime)
  {
    Scope scope("Time Stepping");
    step(test_system, time);
    time += test_system.dt;
    std::cout << "Time: t=" << time << "\t dt=" << test_system.dt << std::endl;
    write_vtk(test_system, time);
  }

  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

  MPI_Finalize();
  LOG::Close();
  Profiler::instance().print();
  return 0;
}
