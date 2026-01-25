#ifndef DENSETREE_H_
#define DENSETREE_H_
#include "zindex.h"
#include <bitset>
#include <cstddef>
#include <utility>
#include <vector>
#define ONES SIZE_MAX
#include <cstdint>
#include <iostream>

namespace DenseTree {

struct DenseTree
{

  uint16_t maxDepth;
  std::vector<size_t> _indices;
  std::vector<uint16_t> _depths;
  std::vector<size_t> _index_cache;
  std::vector<uint16_t> _sizes;

  constexpr void print();
};

constexpr size_t get_sparse_index(DenseTree tree, size_t index)
{
  return tree._index_cache[index];
};

constexpr size_t get_dense_index(DenseTree tree, Zindex index)
{
  size_t idx = 0;
  uint16_t currentDepth = 0;
  for (size_t depth = 0; depth < index.depth; depth++)
  {
    uint16_t subtreeDepth = tree._depths[idx];
    if (subtreeDepth < 1)
    {

      return idx;
    }
    //  Filter out already acounted for depth:
    //  ie. for currentDepth=3 and maxdepth=5 use bitmask (4^3-1)*4^(5-2)=0b11_11_11_00_00
    uint16_t dDepth = index.depth - currentDepth;
    size_t mask = ((size_t(1) << 2 * index.depth) - 1) >> 2 * currentDepth;

    size_t local_index = (index.index & mask) >> 2 * (dDepth - 1);
    idx = tree._indices[idx] + local_index;
    currentDepth += subtreeDepth;
  }

  return idx;
};

void optimize_depths();
void add_layer();
void mark_subtree_unused();

template <typename Operator, typename T, typename... Args>
void mipmap(Operator&& O, DenseTree tree, std::vector<T> _data, Args&&... args)
{
  for (size_t depth = tree.maxDepth - 1; depth > 0; depth--)
  {
    for (size_t local_index = tree._sizes.at(depth - 1); local_index < tree._sizes.at(depth); local_index++)
    {
      if (tree._depths[local_index] > 0)
      {
        size_t data_index = tree._indices[local_index];
        _data[local_index] = std::forward<Operator>(O)({ _data[data_index + 1], _data[data_index + 2], _data[data_index + 3], _data[data_index + 4] }, std::forward<Args>(args)...);
      }
    }
  }
}

template <typename T>
T sum(std::array<T, 4> data)
{
  return data[0] + data[1] + data[2] + data[3];
};
template <typename T>
T mean(std::array<T, 4> data)
{
  return 0.25 * (data[0] + data[1] + data[2] + data[3]);
};
template <typename T>
T max(std::array<T, 4> data)
{
  return 0.25 * (data[0] + data[1] + data[2] + data[3]);
};

// uint16_t has_children(size_t index, size_t height, size_t nx, size_t ny)
//{
//   // size_t x_power = std::bit_width<size_t>(nx - 1);
//   // size_t y_power = std::bit_width<size_t>(ny - 1);
//   //  const size_t size = std::max(x_power, y_power);
//
//   size_t local_cell_size = 1 << height;
//   auto [x, y, depth] = ZorderToIndex(index << height * 2);
//   bool inside = (x + local_cell_size <= nx) && (y + local_cell_size <= ny);
//   bool on_boundary = (x <= nx && nx < x + local_cell_size) || (y <= ny && ny < y + local_cell_size);
//   if ((inside || on_boundary) && (height > 0))
//     return 1;
//   return 0;
// };

template <typename Operator, typename... Args>
DenseTree build_tree(Operator&& O, uint16_t maxDepth, Args&&... args)
{
  std::vector<size_t> _indices;
  std::vector<size_t> _index_cache;
  std::vector<uint16_t> _depth;

  // log2 of the smallest square with size 2^size x 2^size,
  // such that the rectangle with nx x ny fits inside

  std::vector<uint16_t> _sizes;

  // Root node
  _sizes.push_back(0);
  _depth.push_back(1);
  _index_cache.push_back(0);
  _indices.push_back(1);
  _sizes.push_back(1);

  size_t index = 1;
  size_t size = 1;

  // Iterate over depths
  for (uint16_t depth = 0; depth < maxDepth; depth++)
  {
    // Iterate Over All nodes on current depth
    for (size_t local_index = _sizes.at(depth); local_index < _sizes.at(depth + 1); local_index++)
    {
      if (_depth.at(local_index) > 0)
      {
        size_t global_index = _index_cache[local_index];
        for (size_t i = 0; i < 4; i++)
        {
          size_t child_zindex = global_index + (i << 2 * (maxDepth - depth - 1));
          bool has_child = std::forward<Operator>(O)(child_zindex, maxDepth, depth + 1, std::forward<Args>(args)...);
          index += 4 * has_child;
          _indices.push_back(index);
          _depth.push_back(has_child);
          _index_cache.push_back(child_zindex);
          size++;
        }
      }
    }
    _sizes.push_back(size);
  }
  return DenseTree { maxDepth, _indices, _depth, _index_cache, _sizes };
};

template <typename Operator, typename... Args>
void broadcast_subtree(Operator&& O, size_t index, uint8_t depth, DenseTree tree, Args&&... args)
{

  std::forward<Operator>(O)(index, depth, std::forward<Args>(args)...);
  if (tree._depths.at(index) > 0 && depth > 0)
  {
    for (int i = 0; i < 4; i++)
    {
      broadcast_subtree(std::forward<Operator>(O), tree._indices.at(index) + i, depth - 1, tree, std::forward<Args>(args)...);
    }
  }
};
template <typename Operator, typename... Args>
void broadcast_depth_first(Operator&& O, DenseTree tree, uint8_t maxDepth, Args&&... args)
{
  broadcast_subtree(std::forward<Operator>(O), 0, maxDepth, tree, std::forward<Args>(args)...);
};

template <typename Operator, typename... Args>
void broadcast_breath_first(Operator&& O, const DenseTree& tree, uint8_t iter_depth, Args&&... args)
{

  for (uint8_t depth = 0; depth <= iter_depth; depth++)
  {

    for (size_t local_index = tree._sizes.at(depth); local_index < tree._sizes.at(depth + 1); local_index++)
    {
      std::forward<Operator>(O)(local_index, depth, std::forward<Args>(args)...);
    }
  }
};

constexpr void print_node(size_t index, uint8_t depth, const DenseTree& tree)
{
  auto [x, y, d] = ZorderToIndex(tree._index_cache.at(index));

  std::cout << "N" << index << " [label=\"x:" << x << "\ny:" << y << "\nd:" << static_cast<size_t>(tree._depths.at(index)) << "\nh:" << (1 << (tree.maxDepth - depth)) << "\"]" << ";" << std::endl;
  if (tree._depths.at(index) > 0)
  {
    for (size_t i = 0; i < 4; i++)
    {
      size_t local_index = tree._indices.at(index) + i;
      std::cout << "N" << index << "->" << "N" << local_index << ";" << std::endl;
    }
  }
};

constexpr void DenseTree::print()
{

  // std::cout << "Size:\t";
  // for (auto s : _sizes)
  //{
  //   std::cout << static_cast<size_t>(s) << "\t";
  // }
  // std::cout << std::endl;

  // std::cout << "Idx:\t";
  // for (auto idx : _indices)
  //{
  //   std::cout << idx << "\t";
  // }
  // std::cout << std::endl;
  // std::cout << "Cache:\t";
  // for (auto s : _index_cache)
  //{
  //   std::cout << static_cast<size_t>(s) << "\t";
  // }
  // std::cout << std::endl;
  // std::cout << "Depth:\t";
  // for (auto d : _depths)
  //{
  //   std::cout << static_cast<size_t>(d) << "\t";
  // }
  // std::cout << std::endl;
  /*
   * * Prints Tree Structure to dot syntax
   */
  std::cout << "digraph d{\n label=\"visualization of the compressed quadtree of depth: " << static_cast<size_t>(maxDepth) << "\"" << std::endl;
  broadcast_breath_first(print_node, *this, maxDepth, *this);
  std::cout << "}" << std::endl;
};
};

#endif // DENSETREE_H_
