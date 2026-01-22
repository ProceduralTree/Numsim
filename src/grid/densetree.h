#ifndef DENSETREE_H_
#define DENSETREE_H_
#include <algorithm>
#include <bit>
#include <vector>
#define ONES 0xFFFFFFFFF
#include <cstddef>
#include <cstdint>
#include <grid/quadtree.h>

inline uint16_t compact1by1(uint32_t n)
{
  n &= 0x55555555;
  n = (n ^ (n >> 1)) & 0x33333333;
  n = (n ^ (n >> 2)) & 0x0F0F0F0F;
  n = (n ^ (n >> 4)) & 0x00FF00FF;
  n = (n ^ (n >> 8)) & 0x0000FFFF;
  return static_cast<uint16_t>(n);
}

Index ZorderToIndex(uint32_t index)
{
  uint16_t x = compact1by1(index >> 0); // even bits
  uint16_t y = compact1by1(index >> 1); // odd bits
  return Index { x, y };
}

namespace DenseTree {
struct Zindex
{
  size_t index;
  uint8_t depth;
  Zindex(Index I)
    : index(IndexToZOrder(I.x, I.y))
    , depth(I.depth) { };
};

struct DenseTree
{

  size_t _indices[];
  uint8_t _depths[];
  size_t _index_cache[];

  uint8_t _sizes[];
  uint8_t maxDepth;
};

Zindex get_sparse_index(DenseTree tree, size_t index)
{
  return DenseTree._index_cache[index];
};

size_t get_dense_index(DenseTree tree, Zindex index)
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

void optimize_depths();
void add_layer();
void mark_subtree_unused();

template <typename T>
void mipmap(DenseTree tree, std::vector<T> _data)
{
  for (size_t depth = tree.maxDepth - 1; depth > 0; depth--)
  {
    for (size_t local_index = tree._sizes.at(depth - 1); local_index < tree._sizes.at(depth); local_index++)
    {
      if (tree._depths[local_index] > 0)
      {
        size_t data_index = tree._indices[local_index];
        _data[local_index] = 0;
        for (int i = 0; i < 4; i++)
          _data[local_index] += _data[data_index + i];
      }
    }
  }
}

bool has_children(size_t index, size_t height, size_t nx, size_t ny)
{
  size_t local_cell_size = 1 << height;
  auto [x, y] = ZorderToIndex(index);
  if ((x <= nx && nx <= x + local_cell_size) || (y <= ny && ny <= y + local_cell_size))
    return false;
};

void build_from_rectangle(size_t nx, size_t ny)
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

  _sizes.emplace_back(0);
  size_t index = 0;
  // Iterate over depths
  for (uint8_t depth = 0; depth < size; depth++)
  {
    _sizes.emplace_back(index + 1);
    // Iterate Over All nodes on current depth
    for (size_t local_index = _sizes.at(depth - 1); local_index < _sizes.at(depth); local_index++)
    {
      size_t global_index = _index_cache.at(local_index);
      if (has_children(global_index, size - depth, , nx, ny))
      {
        _depth.at(local_index) = 1;
        for (size_t i = 0; i < 4; i++)
        {
          index++;
          _indices.emplace_back(index);
          _depth.emplace_back(0);
          _index_cache.emplace_back(global_index + i);
        }
      }

      // Get maximum subdepth
    }
  }
};

};

#endif // DENSETREE_H_
