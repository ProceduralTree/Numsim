#ifndef UTIL_H_
#define UTIL_H_

#include "grid/densetree.h"
#include <cstddef>
#include <cstdint>
constexpr uint16_t is_desired_depth(size_t index, uint16_t depth, uint16_t maxDepth, uint16_t desired, const DenseTree::DenseTree& old_tree)
{
  if (depth <= desired && depth <= maxDepth)
  {
    return 1;
  }
  // size_t old_index = DenseTree::get_dense_index(old_tree, { index, maxDepth });
  // auto [_, old_depth] = old_tree._index_cache[old_index];
  // if (depth <= old_depth)
  //{
  //   return 1;
  // }
  return 0;
};

#endif // UTIL_H_
