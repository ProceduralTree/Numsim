#ifndef BROADCAST_H_
#define BROADCAST_H_

#include "../pde/system.h"
#include "grid.h"
#include <cstdint>
#include <functional>

void broadcast(std::function<void(PDESystem &, uint16_t, uint16_t)> Operator,
               PDESystem system) {
    // #pragma omp parallel for
    for (uint16_t j = 1; j < system.size_y - 1; j++) {
        for (uint16_t i = 1; i < system.size_x - 1; i++) {
            Operator(system, i, j);
        }
    }
}
#endif // BROADCAST_H_
