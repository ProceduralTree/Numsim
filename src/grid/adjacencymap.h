#ifndef ADJACENCYMAP_H_
#define ADJACENCYMAP_H_

#include "grid/boundary.h"
#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include "grid/zindex.h"
#include "utils/index.h"
#include <cstddef>
#include <cstdint>
constexpr void set_neighbour(size_t index, uint16_t depth, SparseGrid2D<size_t>& grid, Offset o)
{
  Index I = ZorderToIndex(grid.tree._index_cache[index]);
  size_t local_cell_size = 1ULL << (grid.tree.maxDepth - depth);
  Zindex neighbour = IndexToZOrder(I + local_cell_size * o);
  size_t neighbour_index = DenseTree::get_dense_index(grid.tree, neighbour);
  grid[index] = neighbour_index;
  // assert(neighbour.index == grid.tree._index_cache[neighbour_index].index);
}
struct AdjMap
{
  SparseGrid2D<size_t> _top;
  SparseGrid2D<size_t> _bottom;
  SparseGrid2D<size_t> _left;
  SparseGrid2D<size_t> _right;
  AdjMap(const BoundaryFlags& flags)
    : _top(flags.tree)
    , _bottom(flags.tree)
    , _left(flags.tree)
    , _right(flags.tree)
  {
    DenseTree::broadcast_breath_first(set_neighbour, flags.tree, flags.tree.maxDepth, _top, Iy);
    DenseTree::broadcast_breath_first(set_neighbour, flags.tree, flags.tree.maxDepth, _bottom, -Iy);
    DenseTree::broadcast_breath_first(set_neighbour, flags.tree, flags.tree.maxDepth, _left, -Ix);
    DenseTree::broadcast_breath_first(set_neighbour, flags.tree, flags.tree.maxDepth, _right, Ix);
  };
};
template <Offset offset, Sign s>
constexpr size_t get_index(size_t index, const AdjMap& map)
{
  if constexpr (offset == Ix && s == Sign::Plus)
    return map._right[index];
  if constexpr (offset == Ix && s == Sign::Minus)
    return map._left[index];
  if constexpr (offset == Iy && s == Sign::Plus)
    return map._top[index];
  if constexpr (offset == Iy && s == Sign::Minus)
    return map._bottom[index];
};

#endif // ADJACENCYMAP_H_
