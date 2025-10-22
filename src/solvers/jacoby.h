#ifndef JACOBY_H_
#define JACOBY_H_

#include "grid.h"
#include "system.h"

void jacoby(Grid2D &x, const Grid2D &b);
void gaus_seidel_step(PDESystem system);

#endif // JACOBY_H_
