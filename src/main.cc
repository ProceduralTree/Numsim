#include "utils/Logger.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <csignal>
#include <grid/grid.h>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <utils/partitioning.h>
#include <utils/profiler.h>

void signalInt(int sig)
{
  DebugF("Interrupt from: {}", sig);
  Profiler::PrintStack();
  Profiler::Close();
  exit(sig);
}

auto main(int argc, char* argv[]) -> int
{
  signal(SIGINT, signalInt);

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  LOG::Init(LOG::LoggerType::STDOUT);
  Profiler::Init(Profiler::Type::ACCUMULATE);
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
  MPIInfo mpiInfo = MPIInfo();
  setMPIInfo(mpiInfo, Settings::get(), rank, size);
  PDESystem system = PDESystem(Settings::get(), mpiInfo);

  double time = 0;
  ProfilePush("main");
  while (time < system.settings.endTime)
  {
    ProfileCount();
    step(system, time);
    time += system.dt;
    step(system, time);
    time += system.dt;
    for (int i = 0; i < mpiInfo.size; i++)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpiInfo.rank == i)
      {
        std::cout << "Hello from Rank " << rank << " of " << size << std::endl;
        std::cout << "\rTime: t=" << time << "\t dt=" << system.dt << std::endl;
        std::cout << "\rPressure: t=" << system.p << std::endl;
        std::cout << "\r V: " << system.v << "\t U:" << system.u << std::endl;
      }
    }
    //  break;
    // write_vtk(system, time);
  }
  std::cout << std::endl;
  ProfilePop();

  MPI_Finalize();
  LOG::Close();
  Profiler::Close();
  return 0;
}
