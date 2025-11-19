#include "utils/Logger.h"
#include "utils/settings.h"
#include <chrono>
#include <grid/grid.h>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <sstream>
#include <vector>

auto main(int argc, char* argv[]) -> int
{
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  LOG::Init(LOG::LoggerType::STDOUT);
  if (argc < 2)
  {
    LOG::Warning("missing file name");
    LOG::Close();
    return -1;
  }

  if (!Settings::loadFromFile(argv[1]))
  {
    LOG::Warning("couldn't parse settings file");
    LOG::Close();
    return -1;
  }
  Settings::get().printSettings();

  PDESystem system = PDESystem(Settings::get());

  std::vector<std::array<double, 2>> timeSeries;
  double time = 0;
  std::ostringstream oss;
  oss << "data_" << Settings::get().maximumDt << ".csv";
  std::string filename = oss.str();
  while (time < system.settings.endTime && system.p.max() < 1e8)
  {
    step(system, time);
    time += system.dt;
    std::cout << "\rTime: t=" << time << "\t dt=" << system.dt << "\t maxp=" << system.p.max() << std::flush;
    write_vtk(system, time);
    // timeSeries.push_back({ time, system.p.max() });
  }
  std::cout << std::endl;
  std::cout << "Final pressure field:" << std::endl;
  std::cout << system.p << std::endl;
  // std::ofstream file(filename);
  // if (!file)
  // {
  //   std::cerr << "Fehler beim Öffnen der Datei!" << std::endl;
  //   return 1;
  // }
  // file << "Time,MaxPressure" << std::endl;
  // for (auto v : timeSeries)
  // {
  //   file << v[0] << "," << v[1] << "\n";
  // }
  // file << std::flush;
  // file.close();
  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

  MPI_Finalize();
  LOG::Close();
  return 0;
}
