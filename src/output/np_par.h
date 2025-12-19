
#pragma once
#include <utils/partitioning.h>

struct Grid2D;
struct Range;
struct PDESystem;
namespace np_par {

void init(const PDESystem& system);
void writeNP(const PDESystem& system, double dt);
}
