#ifndef DENSETREE_H_
#define DENSETREE_H_
#include <algorithm>
#include <bit>
#include <vector>
#define ONES 0xFFFFFFFFF
#include "quadtree.h"
#include <cstddef>
#include <cstdint>

namespace DenseTree {
struct TreeIndex
{
  size_t index;
  uint8_t depth;
};

struct DenseTree
{
  TreeIndex _data[];
  size_t _indices[];
  uint8_t _depths[];
  TreeIndex _index_cache[];
  uint8_t maxDepth;
};

TreeIndex get_sparse_index(DenseTree tree, size_t index)
{
  return DenseTree._index_cache[index];
};

size_t get_dense_index(DenseTree tree, TreeIndex index)
{
  size_t idx = 0;
  uint8_t currentDepth = 0;
  for (size_t depth = 0; depth < index.depth; depth++)
  {
    uint8_t subtreeDepth = tree._depths[idx];
    uint8_t dDepth = index.depth - subtreeDepth;
    if (dDepth < 1)
      return idx;
    // Filter out already acounted for depth:
    // ie. for currentDepth=3 and maxdepth=5 use bitmask (4^3-1)*4^(5-2)=0b11_11_11_00_00
    size_t mask = ONES >> (sizeof(size_t) * 8 - currentDepth);
    size_t local_index = (index.index & mask) >> 2 * dDepth;
    idx = tree._indices[idx + local_index];
    currentDepth += subtreeDepth;
  }

  return idx;
};

bool has_children(size_t index)
{
  return false;
};

void build_tree_from_settings(size_t nx, size_t ny)
{
  std::vector<size_t> _indices;
  std::vector<size_t> _index_cache;
  std::vector<uint8_t> _depth;

  std::vector<uint8_t> _sizes;
  // log2 of the smallest square with size 2^size x 2^size,
  // such that the rectangle with nx x ny fits inside
  size_t x_power = std::bit_width<size_t>(nx - 1);
  size_t y_power = std::bit_width<size_t>(ny - 1);
  const size_t size = std::max(x_power, y_power);

  // log2 of the larges square with size 2^dense_size x 2^dense_size,
  // such that it fits inside the rectangle nx x ny
  size_t dense_size = (1 << (x_power - 1));
  const size_t leftover = std::min(nx - dense_size, ny - dense_size);
  size_t dense_depth = std::bit_width<size_t>(leftover);

  _sizes.emplace_back(0);
  size_t index = 0;
  // Iterate over depths
  for (uint8_t depth = 0; depth < size; depth++)
  {
    _sizes.emplace_back(index + 1);
    // Iterate Over All nodes on current depth
    for (size_t local_index = _sizes.at(depth - 1); local_index < _sizes.at(depth); local_index++)
    {

      if (has_children(index))
      {
        _indices.emplace_back(index);
        _depth.emplace_back(1);
        index += 4;
      } else
      {
        _indices.emplace_back(0);
        _depth.emplace_back(0);
      }

      // Get maximum subdepth
    }
  }
};

};

#endif // DENSETREE_H_
