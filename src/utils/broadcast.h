#ifndef BROADCAST_H_
#define BROADCAST_H_

#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "grid/zindex.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <cstdint>
#include <grid/grid.h>
#include <pde/system.h>
#include <utils/index.h>

template <typename Operator, typename... Args>
void broadcast_type(Operator&& O, size_t index, uint8_t depth, const BoundaryFlags& flags, uint16_t cell_type, Args&&... args)
{
  bool contains_cell_type = cell_type & flags.flags[index];
  if (contains_cell_type && depth == 0)
  {
    Index I = ZorderToIndex(flags.tree._index_cache[index]);
    std::forward<Operator>(O)(I, std::forward<Args>(args)...);
    return;
  }
  bool has_subtree = flags.tree._depths.at(index) > 0;
  if (contains_cell_type && has_subtree && depth > 0)
  {
    for (int i = 0; i < 4; i++)
    {
      broadcast_type(std::forward<Operator>(O), flags.tree._indices.at(index) + i, depth - 1, flags, cell_type, std::forward<Args>(args)...);
    }
  }
}

template <typename Operator, typename... Args>
void broadcast(Operator&& O, const BoundaryFlags& flags, uint16_t B, Args&&... args)
{
  ProfileScope("Depth First Broadcast");
  broadcast_type(std::forward<Operator>(O), 0, flags.tree.maxDepth, flags, B, std::forward<Args>(args)...);
};

template <typename Operator, typename... Args>
void broadcast_boundary(Operator&& O, const BoundaryFlags& flags, uint16_t B, Args&&... args)
{
  broadcast(std::forward<Operator>(O), flags, static_cast<uint16_t>(BoundaryType::TOP) & B, -Iy, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), flags, static_cast<uint16_t>(BoundaryType::BOTTOM) & B, Iy, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), flags, static_cast<uint16_t>(BoundaryType::LEFT) & B, Ix, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), flags, static_cast<uint16_t>(BoundaryType::RIGHT) & B, -Ix, std::forward<Args>(args)...);
};

inline void copy(Index I, Offset O, const SparseGrid2D<double>& from, SparseGrid2D<double>& to) { to[I] = from[I]; };
// inline void copy(Index I, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };

inline void set(Index I, Offset O, SparseGrid2D<double>& array, double value) { array[I] = value; };
#endif // BROADCAST_H_
