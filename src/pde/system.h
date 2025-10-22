#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "grid.h"
#include "indexing.h"

struct PDESystem {
    Grid2D p;
    Grid2D u;
    Grid2D v;
    Grid2D F;
    Grid2D G;
    Grid2D b;
    uint16_t size_x;
    uint16_t size_y;
};

#endif // SYSTEM_H_
