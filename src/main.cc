#include "grid/boundary.h"
#include "grid/densetree.h"
#include "grid/util.h"
#include "output/vtk_tree.h"
#include "pde/pressuresolvers.h"
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
#include <pde/system.h>
#include <sstream>
#include <utils/profiler.h>

void signalInt(int sig)
{
  DebugF("Interrupt from: {}", sig);
  Profiler::Close();
  exit(sig);
}

auto main(int argc, char* argv[]) -> int
{
  signal(SIGINT, signalInt);
  signal(SIGTERM, signalInt);

  LOG::Init(LOG::LoggerType::FILE);
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
  // Settings::get().printSettings();
  auto r = Range { Index { 1, 1, 0 }, Index { static_cast<uint16_t>(Settings::get().nCells[0] + 1), static_cast<uint16_t>(Settings::get().nCells[1] + 1), 0 } };
  auto t0 = DenseTree::from_range(r);

  auto t = dilate(t0);
  auto flags0 = BoundaryFlags(t, r);
  auto tree = dilate(flags0);
  // auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth - 1, t);
  //  updated_tree.print();
  auto flags = BoundaryFlags(tree, r);
  PDESystem system = PDESystem(Settings::get(), flags);
  auto solver = Jacoby(tree);

  double time = 0;

  double next_written_time = 0;
  std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
  // std::chrono::system_clock::time_point last_time = std::chrono::system_clock::now();

  ProfileScope("main");
  while (time < system.settings.endTime)
  {
    // if (system.dt < 1e-16)
    //{
    //   std::cerr << "To Small TimeStep" << std::endl;
    //   abort();
    // }
    step(system, solver, time);
    time += system.dt;

    if (time > next_written_time)
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

      fflush(stdout);
      next_written_time += 1;
      auto data_set = init(tree, false);
      write_field("Boundary Data", system.boundary.tree, system.boundary.flags._data, data_set);
      write_field("U raw data", system.boundary.tree, system.u._data, data_set);
      write_field("V raw data", system.boundary.tree, system.v._data, data_set);
      write_field("P raw data", system.boundary.tree, system.p._data, data_set);
      write_field("F raw data", system.boundary.tree, system.F._data, data_set);
      write_field("G raw data", system.boundary.tree, system.G._data, data_set);
      write_field("RHS raw data", system.boundary.tree, system.rhs._data, data_set);
      write_field("Residual", system.boundary.tree, solver.residual._data, data_set);
      save_dataset(data_set);
    }
    // auto data_set = init(tree, true);
    // write("Solve", system, data_set);
    // save_dataset(data_set);
  }
  std::cout << std::endl;

  LOG::Close();
  Profiler::Close();
  return 0;
}
