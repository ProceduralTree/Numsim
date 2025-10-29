#include "settings.h"

#include <cassert>
#include <cstdlib>
#include <sstream>

#define LOG(...) 0

static std::string testString = ""
                                "# Settings file for numsim program\n"
                                "# Run ./numsim lid_driven_cavity.txt\n"
                                "\n"
                                "# Problem description\n"
                                "physicalSizeX = 2.0   # physical size of the domain\n"
                                "physicalSizeY = 2.0\n"
                                "endTime = 10.0        # duration of the simulation\n"
                                "re = 1000             # Reynolds number\n"
                                "gX = 0.0              # external forces, set to (gX,gY) = (0,-9.81) to account for gravity\n"
                                "gY = 0.0\n"
                                "\n"
                                "# Dirichlet boundary conditions\n"
                                "dirichletBottomX = 0\n"
                                "dirichletBottomY = 0\n"
                                "dirichletTopX    = 1\n"
                                "dirichletTopY    = 0\n"
                                "dirichletLeftX   = 0\n"
                                "dirichletLeftY   = 0\n"
                                "dirichletRightX  = 0\n"
                                "dirichletRightY  = 0\n"
                                "\n"
                                "# Discretization parameters\n"
                                "nCellsX = 20          # number of cells in x and y direction\n"
                                "nCellsY = 20\n"
                                "useDonorCell = true   # if donor cell discretization should be used, possible values: true false\n"
                                "alpha = 0.5           # factor for donor-cell scheme, 0 is equivalent to central differences\n"
                                "tau = 0.5             # safety factor for time step width\n"
                                "maximumDt = 0.1       # maximum values for time step width\n"
                                "\n"
                                "# Solver parameters\n"
                                "pressureSolver = SOR  # which pressure solver to use, possible values: GaussSeidel SOR CG\n"
                                "omega = 1.6           # overrelaxation factor, only for SOR solver\n"
                                "epsilon = 1e-5        # tolerance for 2-norm of residual\n"
                                "maximumNumberOfIterations = 1e4    # maximum number of iterations in the solver";

bool compareToSecond(const char* a, const char* b)
{
  while (*a == *b)
  {
    if (*a == '\0')
      return *b == '\0';
    if (*b == '\0')
      return true;
    a++;
    b++;
  }
  return false;
}
void parseSettings(const std::string& input, Settings* settings)
{
  std::stringstream stream(input);
  uint32_t lineNumber = 0;
  for (std::string line; std::getline(stream, line); lineNumber++)
  {
    if (line.starts_with('#'))
    {
      continue;
    }
    bool validLine = line.contains('=');
    if (!validLine)
    {
      // TODO: use logger or discard error
      LOG("invalid Line found at : %d", lineNumber);
      continue;
    }
    size_t eqPos = line.find_first_of("=") + 1;
    size_t maxLength = line.length();
    size_t keyoffset = line.find_first_not_of(" \t");
    size_t keylength = line.find_first_of("= \t") - keyoffset;
    assert(keyoffset + keylength < maxLength && "error reading key");
    const char* key = line.c_str() + keyoffset;
    size_t valueoffset = line.find_first_not_of(" \t", eqPos);
    assert(valueoffset < maxLength && "error reading value");
    std::string value = line.substr(valueoffset, maxLength - valueoffset);
    if (compareToSecond(key, "nCells"))
    {
      if (compareToSecond(key + sizeof("nCells") - 1, "X"))
        settings->nCells[0] = atoi(value.c_str());
      else if (compareToSecond(key + sizeof("nCells") - 1, "Y"))
        settings->nCells[1] = atoi(value.c_str());
      else
        validLine = false;
    } else if (compareToSecond(key, "physicalSize"))
    {
      if (compareToSecond(key + sizeof("physicalSize") - 1, "X"))
        settings->physicalSize[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("physicalSize") - 1, "Y"))
        settings->physicalSize[1] = atof(value.c_str());
      else
        validLine = false;
    } else if (compareToSecond(key, "re"))
      settings->re = atof(value.c_str());
    else if (compareToSecond(key, "endTime"))
      settings->endTime = atof(value.c_str());
    else if (compareToSecond(key, "tau"))
      settings->tau = atof(value.c_str());
    else if (compareToSecond(key, "maximumDt"))
      settings->maximumDt = atof(value.c_str());
    else if (compareToSecond(key, "g"))
    {
      if (compareToSecond(key + sizeof("g") - 1, "X"))
        settings->g[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("g") - 1, "Y"))
        settings->g[1] = atof(value.c_str());
      else
        validLine = false;
    } else if (compareToSecond(key, "useDonorCell"))
      settings->useDonorCell = value.starts_with("true");
    else if (compareToSecond(key, "alpha"))
      settings->alpha = atof(value.c_str());
    else if (compareToSecond(key, "dirichletBc"))
    {
      if (compareToSecond(key + sizeof("dirichlet") - 1, "BottomX"))
        settings->dirichletBcBottom[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "BottomY"))
        settings->dirichletBcBottom[1] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "TopX"))
        settings->dirichletBcTop[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "TopY"))
        settings->dirichletBcTop[1] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "LeftX"))
        settings->dirichletBcLeft[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "LeftY"))
        settings->dirichletBcLeft[1] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "RightX"))
        settings->dirichletBcRight[0] = atof(value.c_str());
      else if (compareToSecond(key + sizeof("dirichlet") - 1, "RightY"))
        settings->dirichletBcRight[1] = atof(value.c_str());
      else
        validLine = false;
    } else if (compareToSecond(key, "pressureSolver"))
    {
      if (value.starts_with("SOR"))
        settings->pressureSolver = Settings::PressureSolver::SOR;
      else if (value.starts_with("GaussSeidel"))
        settings->pressureSolver = Settings::PressureSolver::GaussSeidel;
      else
        validLine = false;
    } else if (compareToSecond(key, "omega"))
      settings->omega = atof(value.c_str());
    else if (compareToSecond(key, "epsilon"))
      settings->epsilon = atof(value.c_str());
    else if (compareToSecond(key, "maximumNumberOfIterations"))
      settings->maximumNumberOfIterations = atof(value.c_str());
    else
      validLine = false;
    // once again check for invalid line
    if (!validLine)
      LOG("invalid line found at : %d", invalidLine);
  }
}

void Settings::loadFromFile(std::string filename)
{
}

void Settings::printSettings()
{
}
