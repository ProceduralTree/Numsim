#pragma once
#include <utils/partitioning.h>

struct Grid2D;
struct Range;
struct PDESystem;
namespace vtk_par {

void init(const PDESystem& system);
void writeVTK(const PDESystem& system, double dt);
}
