#include "utils/Logger.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <grid/grid.h>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <vector>

auto main(int argc, char* argv[]) -> int
{
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  LOG::Init(LOG::LoggerType::FILE);
  Profiler::Init();
  if (argc < 2)
  {
    LOG::Warning("missing file name");
    LOG::Close();
    Profiler::Close();
    return -1;
  }

  if (!Settings::loadFromFile(argv[1]))
  {
    LOG::Warning("couldn't parse settings file");
    LOG::Close();
    Profiler::Close();
    return -1;
  }
  Settings::get().printSettings();

  PDESystem system = PDESystem(Settings::get());

  double time = 0;
  ProfilePush("main");
  while (time < system.settings.endTime)
  {
    ProfileCount();
    step(system, time);
    time += system.dt;
    std::cout << "\rTime: t=" << time << "\t dt=" << system.dt << std::flush;
    write_vtk(system, time);
  }
  std::cout << std::endl;
  ProfilePop();

  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

  MPI_Finalize();
  LOG::Close();
  Profiler::Close();
  return 0;
}
