#include "broadcast.h"
#include "grid.h"
#include "indexing.h"

void broadcast(std::function<void(PDESystem &, unit16_t, uint16_t)> Operator,
               PDESystem system) {
#pragma omp parallel for
    for (uint32_t j = 1; j < system.size_y - 1; j++) {
        for (uint32_t i = 1; i < system.size_x - 1; i++) {
            Operator(system, i, j);
        }
    }
}

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
