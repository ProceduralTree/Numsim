#ifndef RECTANGLE_H_
#define RECTANGLE_H_
/*
 * * Provides Util Functions,
 * to enable the creation of a rectangular Domain on a sparse quadtree grid
 */

#include "grid/zindex.h"
#include "utils/index.h"
#include <cstdint>
#include <format>
#include <iostream>

constexpr uint16_t intersects_range(size_t index, uint16_t depth, Range r, uint16_t maxDepth)
{
  if (depth >= maxDepth)
    return 0;

  size_t local_cell_size = 1 << (maxDepth - depth);
  auto [x, y, _] = ZorderToIndex(index);
  // std::cerr << std::format("Index({},{}) \t, local cell = {} , Range(({},{})->({},{}))", x, y, local_cell_size, r.begin.x, r.begin.y, r.end.x, r.end.y) << std::endl;
  bool on_right_boundary = (x <= r.end.x + local_cell_size && r.end.x < x + 2 * local_cell_size);
  bool on_left_boundary = (x <= r.begin.x && r.begin.x < x + 2 * local_cell_size);
  bool on_bottom_boundary = (y <= r.begin.y && r.begin.y < y + 2 * local_cell_size);
  bool on_top_boundary = (y <= r.end.y + local_cell_size && r.end.y < y + 2 * local_cell_size);
  bool in_range = x <= r.end.x + local_cell_size && y <= r.end.y + local_cell_size;
  // bool in_range = true;
  if (on_top_boundary | on_bottom_boundary | on_left_boundary | on_right_boundary && in_range)
    return 1;
  return 0;
};

constexpr uint16_t intersects_top(size_t index, uint16_t depth, Range r, uint16_t maxDepth)
{
  if (depth >= maxDepth)
    return 0;

  size_t local_cell_size = 1 << (maxDepth - depth);
  auto [x, y, _] = ZorderToIndex(index);
  // std::cerr << std::format("Index({},{}) \t, local cell = {} , Range(({},{})->({},{}))", x, y, local_cell_size, r.begin.x, r.begin.y, r.end.x, r.end.y) << std::endl;
  bool on_right_boundary = (x <= r.end.x && r.end.x < x + local_cell_size);
  bool on_top_boundary = (y <= r.end.y && r.end.y < y + local_cell_size);
  bool in_range = x <= r.end.x && y <= r.end.y;
  // bool in_range = true;
  if (on_top_boundary | on_right_boundary && in_range)
    return 1;
  return 0;
};
#endif // RECTANGLE_H_
