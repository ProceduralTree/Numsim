#include "output/vtk_par.h"
#include "utils/Logger.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <grid/grid.h>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <sstream>
#include <utils/partitioning.h>
#include <utils/profiler.h>

void signalInt(int sig)
{
  DebugF("Interrupt from: {}", sig);
  Profiler::Close();
  MPI_Barrier(MPI_COMM_WORLD);
  exit(sig);
}

auto main(int argc, char* argv[]) -> int
{
  signal(SIGINT, signalInt);
  signal(SIGTERM, signalInt);

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  LOG::Init(LOG::LoggerType::STDOUT);
  Profiler::Init(Profiler::Type::ACCUMULATEPAR);
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
  // Settings::get().printSettings();
  Partitioning::MPIInfo mpiInfo = Partitioning::MPIInfo();
  setMPIInfo(mpiInfo, Settings::get(), rank, size);
  Settings::set().mpi = mpiInfo;
  PDESystem system = PDESystem(Settings::get(), mpiInfo);

  vtk_par::init(system);

  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;
  std::cout << "nX " << mpiInfo.nCells[0] << " nY " << mpiInfo.nCells[1] << std::endl;

#define DebugPrintGrid(name, g) DebugF("Rank{}: Grid " #name " with size:{},{}, begin:{},{}, end:{},{}, range:({},{}),({},{}), globalRange:({},{}),({},{})", mpiInfo.rank, g.size_x, g.size_y, g.begin.x, g.begin.y, g.end.x, g.end.y, g.range.begin.x, g.range.begin.y, g.range.end.x, g.range.end.y, g.globalRange.begin.x, g.globalRange.begin.y, g.globalRange.end.x, g.globalRange.end.y);

  DebugPrintGrid(p, system.p);
  DebugPrintGrid(u, system.u);
  DebugPrintGrid(v, system.v);
  double time = 0;

  double next_written_time = 1;
  std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
  std::chrono::system_clock::time_point last_time = std::chrono::system_clock::now();

  ProfileScope("main");
  while (time < system.settings.endTime)
  {
    // if (system.dt < 1e-16)
    //{
    //   std::cerr << "To Small TimeStep" << std::endl;
    //   abort();
    // }
    step(system, time);
    time += system.dt;
    step(system, time);
    time += system.dt;
    if (time > next_written_time)
    {
      if (mpiInfo.rank == 0)
      {
        std::chrono::system_clock::time_point tmp_time = std::chrono::system_clock::now();
        auto diff = tmp_time - start_time;
        std::stringstream s;
        s << "\r[";
        for (int i = 0; i < Settings::get().endTime; i++)
        {
          s << ((i < time) ? '#' : ' ');
        }
        s << "]";
        s << "\t Time:" << time << "/" << Settings::get().endTime << "s";
        s << "\t Iter/s:" << std::chrono::duration<double>(diff).count() / time;
        s << "\t Wall Time:" << std::chrono::duration<double>(diff).count();
        s << "\n";
        printf("%s", s.str().c_str());

        fflush(stdout);
      }
      vtk_par::writeVTK(system, time);

      next_written_time++;
    }
    // write_vtk(system, time);
  }
  std::cout << std::endl;

  MPI_Finalize();
  LOG::Close();
  Profiler::Close();
  return 0;
}
