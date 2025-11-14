#pragma once

#include <filesystem>

struct Settings
{
  Settings() = default;
  Settings(std::filesystem::path filename) { loadFromFile(filename); }

  int nCells[2]; //< number of cells in x and y direction
  double physicalSize[2]; //< physical size of the domain
  double re = 400; //< reynolds number
  double endTime = 10.0; //< end time of the simulation
  double tau = 0.8; //< safety factor for time step width
  double maximumDt = 0.1; //< maximum time step width

  double g[2]; //< external forces

  bool useDonorCell = false; //< if the donor cell scheme schould be used
  double alpha = 0.4; //< factor for donor-cell scheme

  double dirichletBcBottom[2]; //< prescribed values of u,v at bottom of domain
  double dirichletBcTop[2]; //< prescribed values of u,v at top of domain
  double dirichletBcLeft[2]; //< prescribed values of u,v at left of domain
  double dirichletBcRight[2]; //< prescribed values of u,v at right of domain

  enum PressureSolver
  {
    BlackRed,
    CG,
    SOR,
    Jacoby,
    GaussSeidel
  };
  PressureSolver pressureSolver = SOR; //< which pressure solver to use, "GaussSeidel" or "SOR"
  double omega = 1.8; //< overrelaxation factor
  double epsilon = 1e-4; //< tolerance for the residual in the pressure solver
  int maximumNumberOfIterations = 1e5; //< maximum number of iterations in the solver

  //! parse a text file with settings, each line contains "<parameterName> = <value>"
  static bool loadFromFile(std::filesystem::path filename);

  //! output all settings to console
  void printSettings() const;
  static const Settings& get();

private:
  bool loadFromFileInternal(std::filesystem::path filename);
};
