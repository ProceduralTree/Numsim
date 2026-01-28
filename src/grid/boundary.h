#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "grid/zindex.h"
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <grid/densetree.h>
#include <grid/sparsegrid.h>

enum class BoundaryType : uint16_t
{
  // sides
  U_BOTTOM = 0b0000'0000'0000'0001,
  U_TOP = 0b0000'0000'0000'0010,
  U_LEFT = 0b0000'0000'0000'0100,
  U_RIGHT = 0b0000'0000'0000'1000,
  U_Inside = 0b0000'0000'0001'0000,
  U = U_BOTTOM | U_LEFT | U_RIGHT | U_TOP | U_Inside,
  U_BOUNDARY = U_BOTTOM | U_LEFT | U_RIGHT | U_TOP,
  P_BOTTOM = 0b0000'0000'0010'0000,
  P_TOP = 0b0000'0000'0100'0000,
  P_LEFT = 0b0000'0000'1000'0000,
  P_RIGHT = 0b0000'0001'0000'0000,
  P_Inside = 0b0000'0010'0000'0000,
  P = P_BOTTOM | P_LEFT | P_RIGHT | P_TOP | P_Inside,
  P_BOUNDARY = P_BOTTOM | P_LEFT | P_RIGHT | P_TOP,
  V_BOTTOM = 0b0000'0100'0000'0000,
  V_TOP = 0b0000'1000'0000'0000,
  V_LEFT = 0b0001'0000'0000'0000,
  V_RIGHT = 0b0010'0000'0000'0000,
  V_Inside = 0b0100'0000'0000'0000,
  V = V_BOTTOM | V_LEFT | V_RIGHT | V_TOP | V_Inside,
  V_BOUNDARY = V_BOTTOM | V_LEFT | V_RIGHT | V_TOP,
};

constexpr void set_u_boundary(size_t local_index, uint16_t depth, SparseGrid2D<uint16_t>& flags)
{
  auto I = ZorderToIndex(flags.tree._index_cache[local_index].index);
  // I.depth = flags.tree.maxDepth;
  I.depth = depth;
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_Inside) && !(flags[I + Ix] & static_cast<uint16_t>(BoundaryType::P_RIGHT)))
    flags[local_index] = static_cast<size_t>(BoundaryType::U_Inside) | flags[local_index];

  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_TOP))
    flags[local_index] = static_cast<size_t>(BoundaryType::U_TOP) | flags[local_index];
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_BOTTOM))
    flags[local_index] = static_cast<size_t>(BoundaryType::U_BOTTOM) | flags[local_index];
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_LEFT))
    flags[local_index] = static_cast<size_t>(BoundaryType::U_LEFT) | flags[local_index];

  if (flags[I + Ix] & static_cast<uint16_t>(BoundaryType::P_RIGHT))
    flags[local_index] = static_cast<size_t>(BoundaryType::U_RIGHT) | flags[local_index];
};
constexpr void set_v_boundary(size_t local_index, uint16_t depth, SparseGrid2D<uint16_t>& flags)
{
  auto I = ZorderToIndex(flags.tree._index_cache[local_index].index);
  // I.depth = flags.tree.maxDepth;
  I.depth = depth;
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_Inside) && !(flags[I + Iy] & static_cast<uint16_t>(BoundaryType::P_TOP)))
    flags[local_index] = static_cast<size_t>(BoundaryType::V_Inside) | flags[local_index];

  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_RIGHT))
    flags[local_index] = static_cast<size_t>(BoundaryType::V_RIGHT) | flags[local_index];
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_BOTTOM))
    flags[local_index] = static_cast<size_t>(BoundaryType::V_BOTTOM) | flags[local_index];
  if (flags[local_index] & static_cast<uint16_t>(BoundaryType::P_LEFT))
    flags[local_index] = static_cast<size_t>(BoundaryType::V_LEFT) | flags[local_index];

  if (flags[I + Iy] & static_cast<uint16_t>(BoundaryType::P_TOP))
    flags[local_index] = static_cast<size_t>(BoundaryType::V_TOP) | flags[local_index];
};

constexpr void set_rectangle_p_boundary_type(size_t local_index, uint16_t depth, Range r, SparseGrid2D<uint16_t>& flags)
{
  auto I = ZorderToIndex(flags.tree._index_cache[local_index].index);
  bool in_p_boundary = I >= r.begin && I <= r.end;
  if (in_p_boundary)
  {
    flags[local_index] = static_cast<uint16_t>(BoundaryType::P_Inside);
    return;
  }
  if (I >= r.begin - Ix && I <= r.end)
  {
    flags[local_index] = static_cast<uint16_t>(BoundaryType::P_LEFT);
    return;
  }
  if (I >= r.begin && I <= r.end + Iy)
  {
    flags[local_index] = static_cast<uint16_t>(BoundaryType::P_TOP);
    return;
  }
  if (I >= r.begin - Iy && I <= r.end)
  {
    flags[local_index] = static_cast<uint16_t>(BoundaryType::P_BOTTOM);
    return;
  }
  if (I >= r.begin && I <= r.end + Ix)
  {
    flags[local_index] = static_cast<uint16_t>(BoundaryType::P_RIGHT);
    return;
  }
};

struct BoundaryFlags
{
  const DenseTree::DenseTree& tree;
  SparseGrid2D<uint16_t> flags;
  BoundaryFlags(const DenseTree::DenseTree& tree, Range pressure_range)
    : tree(tree)
    , flags(tree)
  {
    DenseTree::broadcast_breath_first(set_rectangle_p_boundary_type, tree, tree.maxDepth, pressure_range, flags);
    DenseTree::broadcast_breath_first(set_u_boundary, tree, tree.maxDepth, flags);
    DenseTree::broadcast_breath_first(set_v_boundary, tree, tree.maxDepth, flags);
    mipmap(_or<uint16_t>, flags);
  };
};

template <typename Operator, typename... Args>
void broadcast_cell_type(Operator&& O, size_t index, uint8_t depth, const BoundaryFlags& flags, uint16_t cell_type, Args&&... args)
{
  bool contains_cell_type = cell_type & flags.flags[index];
  if (depth == 0 && contains_cell_type)
  {
    std::forward<Operator>(O)(index, depth, std::forward<Args>(args)...);
    return;
  }
  bool has_subtree = flags.tree._depths.at(index) > 0;
  if (contains_cell_type && has_subtree && depth > 0)
  {
    for (int i = 0; i < 4; i++)
    {
      broadcast_cell_type(std::forward<Operator>(O), flags.tree._indices.at(index) + i, depth - 1, flags, cell_type, std::forward<Args>(args)...);
    }
  }
}

template <typename Operator, typename... Args>
void tree_broadcast(Operator&& O, const BoundaryFlags& flags, uint16_t B, Args&&... args)
{
  broadcast_cell_type(std::forward<Operator>(O), 0, flags.tree.maxDepth, flags, B, std::forward<Args>(args)...);
};
#endif // BOUNDARY_H_
