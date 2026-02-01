#ifndef DENSETREE_H_
#define DENSETREE_H_
#include "grid/rectangle.h"
#include "utils/profiler.h"
#include "zindex.h"
#include <cassert>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>
#define ONES SIZE_MAX
#include <cstdint>
#include <iostream>

namespace DenseTree {

struct DenseTree;
template <typename Operator, typename... Args>
DenseTree build_tree(Operator&& has_children, uint16_t maxDepth, Args&&... args);

struct DenseTree
{

  uint16_t maxDepth;
  std::vector<size_t> _indices;
  std::vector<uint16_t> _depths;
  std::vector<Zindex> _index_cache;
  std::vector<size_t> _sizes;
  // Delete copy constructor and copy assignment
  DenseTree(const DenseTree&) = delete;
  DenseTree& operator=(const DenseTree&) = delete;

  // move constructor and move assignment
  DenseTree(DenseTree&&) = default;
  DenseTree& operator=(DenseTree&&) = default;
  // default constructor
  DenseTree() = default;

  constexpr void print();
};

constexpr DenseTree from_range(Range r)
{
  size_t x_power = std::bit_width<size_t>(r.end.x - 1);
  size_t y_power = std::bit_width<size_t>(r.end.y - 1);
  const uint16_t maxDepth = std::max(x_power, y_power);
  return build_tree(intersects_range, maxDepth, r, maxDepth);
};

constexpr Zindex get_sparse_index(DenseTree tree, size_t index)
{
  return tree._index_cache[index];
};

constexpr size_t get_dense_index(const DenseTree& tree, Zindex index)
{
  ProfileScope("Index Lookup");
  size_t idx = 0;
  uint16_t currentDepth = 0;
  for (size_t depth = 0; depth <= index.depth; depth++)
  {
    uint16_t subtreeDepth = tree._depths[idx];
    if (subtreeDepth == 0 || currentDepth == index.depth || currentDepth == tree.maxDepth)
      return idx;
    //  Filter out already acounted for depth:
    //  ie. for currentDepth=3 and maxdepth=5 use bitmask (4^3-1)*4^(5-2)=0b11_11_11_00_00
    uint16_t dDepth = tree.maxDepth - currentDepth;
    size_t mask = ((size_t(1) << (2 * tree.maxDepth)) - 1) >> (2 * currentDepth);

    size_t local_index = (index.index & mask) >> 2 * (dDepth - 1);
    size_t new_index = tree._indices.at(idx) + local_index;
    // idx = (new_index < tree._sizes[currentDepth]) ? new_index : idx;
    idx = new_index;
    currentDepth += subtreeDepth;
  }

  assert(false);

  return idx;
};

void optimize_depths();
void add_layer();
void mark_subtree_unused();

template <typename Operator, typename... Args>
DenseTree build_tree(Operator&& has_children, uint16_t maxDepth, Args&&... args)
{
  DenseTree tree;
  tree.maxDepth = maxDepth;
  // log2 of the smallest square with size 2^size x 2^size,
  // such that the rectangle with nx x ny fits inside

  // Root node
  tree._sizes.push_back(0);
  tree._depths.push_back(1);
  tree._index_cache.push_back({ 0, 0 });
  tree._indices.push_back(1);
  tree._sizes.push_back(1);

  size_t index = 1;
  size_t size = 1;

  // Iterate over depths
  for (uint16_t depth = 0; depth < maxDepth; depth++)
  {
    // Iterate Over All nodes on current depth
    for (size_t local_index = tree._sizes.at(depth); local_index < tree._sizes.at(depth + 1); local_index++)
    {
      if (tree._depths.at(local_index) > 0)
      {
        auto [global_index, _] = tree._index_cache[local_index];
        for (size_t i = 0; i < 4; i++)
        {

          size_t child_zindex = global_index + (i << (2 * static_cast<size_t>(maxDepth - depth - 1)));
          uint16_t has_child = std::forward<Operator>(has_children)(child_zindex, depth + 1, std::forward<Args>(args)...);
          has_child = has_child && ((depth + 1) <= maxDepth);
          index += 4 * has_child;
          tree._indices.push_back(has_child ? index : SIZE_MAX);
          tree._depths.push_back(has_child);
          tree._index_cache.push_back({ child_zindex, static_cast<uint16_t>(1 + depth) });
          size++;
        }
      }
    }
    tree._sizes.push_back(size);
  }
  return tree;
};

template <typename Operator, typename... Args>
void broadcast_subtree(Operator&& O, size_t index, uint8_t depth, const DenseTree& tree, Args&&... args)
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
void broadcast_depth_first(Operator&& O, const DenseTree& tree, uint8_t maxDepth, Args&&... args)
{
  broadcast_subtree(std::forward<Operator>(O), 0, maxDepth, tree, std::forward<Args>(args)...);
};

template <typename Operator, typename... Args>
constexpr void broadcast_level(Operator&& O, const DenseTree& tree, uint16_t depth, Args&&... args)
{

  for (size_t local_index = tree._sizes.at(depth); local_index < tree._sizes.at(depth + 1); local_index++)
  {
    std::forward<Operator>(O)(local_index, depth, std::forward<Args>(args)...);
  }
}

template <typename Operator, typename... Args>
void broadcast_breath_first(Operator&& O, const DenseTree& tree, uint16_t iter_depth, Args&&... args)
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
  auto [x, y, d] = ZorderToIndex(tree._index_cache.at(index).index);

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
  /*
   * * Prints Tree Structure to dot syntax
   */
  std::cout << "digraph d{\n label=\"visualization of the compressed quadtree of depth: " << static_cast<size_t>(maxDepth) << "\"" << std::endl;
  broadcast_breath_first(print_node, *this, maxDepth, *this);
  std::cout << "}" << std::endl;
};
};
#endif // DENSETREE_H_
