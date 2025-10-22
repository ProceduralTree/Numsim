#include "broadcast.h"
#include "grid.h"
#include "indexing.h"

template <typename... Args>
void broadcast<Indexing::Cartesian>(
    std::function<void(const Grid2D &, Grid2D &)> Operator, Args... args) {

#pragma omp parallel for
    for (uint32_t j = 1; j < in.size_y - 1; j++) {
        for (uint32_t i = 1; i < in.size_x - 1; i++) {
            Operator(i, j, std::forward<Args>(args)...);
        }
    }
}

template <typename, typename... Args>
void broadcast<Indexing::ZOrder>(
    std::function<void(const Grid2D<Indexing::ZOrder> &,
                       Grid2D<Indexing::ZOrder> &)>
        Operator,
    Args... args) {

#pragma omp parallel for
    for (uint32_t z = 0; z < in.elements(); z++) {
        Operator(z, std::forward<Args>(args)...);
    }
}
