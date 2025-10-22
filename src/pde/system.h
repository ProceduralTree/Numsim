#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "grid.h"
#include "indexing.h"

template <Indexing I> struct PDESystem {
    Grid2D<I> p;
    Grid2D<I> u;
    Grid2D<I> v;
    Grid2D<I> F;
    Grid2D<I> G;
    uint_16_t size_x;
    uint_16_t size_y;
};

#endif // SYSTEM_H_
