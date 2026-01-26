#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "grid/zindex.h"
#include <cstddef>
#include <cstdint>
#include <grid/densetree.h>

enum class BoundaryType : uint16_t
{
  // sides
  U_BOTTOM = 0b0000'0000'0000'0001,
  U_TOP = 0b0000'0000'0000'0010,
  U_LEFT = 0b0000'0000'0000'0100,
  U_RIGHT = 0b0000'0000'0000'1000,
  U_Inside = 0b0000'0000'0001'0000,
  P_BOTTOM = 0b0000'0000'0010'0000,
  P_TOP = 0b0000'0000'0100'0000,
  P_LEFT = 0b0000'0000'1000'0000,
  P_RIGHT = 0b0000'0001'0000'0000,
  P_Inside = 0b0000'0010'0000'0000,
  V_BOTTOM = 0b0000'0100'0000'0000,
  V_TOP = 0b0000'1000'0000'0000,
  V_LEFT = 0b0001'0000'0000'0000,
  V_RIGHT = 0b0010'0000'0000'0000,
  V_Inside = 0b0100'0000'0000'0000,
};

constexpr void set_u_boundary(size_t local_index, uint16_t depth, std::vector<uint16_t>& flags, const DenseTree::DenseTree& tree)
{
  if (flags.at(local_index) & BoundaryType::P_Inside)
  {
    flags.at(local_index) |= static_cast<size_t>(BoundaryType::U_Inside);
  }
};

constexpr void set_v_boundary(size_t local_index, uint16_t depth, std::vector<uint16_t>& flags, const DenseTree::DenseTree& tree) {

};

constexpr void set_rectangle_p_boundary_type(size_t local_index, uint16_t depth, Range r, std::vector<uint16_t>& flags, const DenseTree::DenseTree& tree)
{
  auto I = ZorderToIndex(tree._index_cache[local_index].index);
  bool in_p_boundary = I >= r.begin && I <= r.end;
  if (in_p_boundary)
  {
    flags.at(local_index) = static_cast<uint16_t>(BoundaryType::P_Inside);
    return;
  }
  if (I >= r.begin - Ix && I <= r.end)
  {
    flags.at(local_index) = static_cast<uint16_t>(BoundaryType::P_LEFT);
    return;
  }
  if (I >= r.begin && I <= r.end + Iy)
  {
    flags.at(local_index) = static_cast<uint16_t>(BoundaryType::P_TOP);
    return;
  }
  if (I >= r.begin - Iy && I <= r.end)
  {
    flags.at(local_index) = static_cast<uint16_t>(BoundaryType::P_BOTTOM);
    return;
  }
  if (I >= r.begin && I <= r.end + Ix)
  {
    flags.at(local_index) = static_cast<uint16_t>(BoundaryType::P_RIGHT);
    return;
  }
};

struct BoundaryFlags
{
  const DenseTree::DenseTree& tree;
  std::vector<uint16_t> boundary_flags;
  BoundaryFlags(const DenseTree::DenseTree& tree, Range pressure_range)
    : tree(tree)
    , boundary_flags(tree._sizes.at(tree.maxDepth + 1))
  {
    std::cerr << "Begin boundary build" << std::endl;
    DenseTree::broadcast_breath_first(set_rectangle_p_boundary_type, tree, tree.maxDepth, pressure_range, boundary_flags, tree);
    DenseTree::broadcast_breath_first(set_u_boundary, tree, tree.maxDepth, boundary_flags, tree);
    DenseTree::broadcast_breath_first(set_v_boundary, tree, tree.maxDepth, boundary_flags, tree);
  };
};

#endif // BOUNDARY_H_
