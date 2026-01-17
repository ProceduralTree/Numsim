#ifndef QUADTREE_H_
#define QUADTREE_H_
#include "utils/index.h"
#include <cstdint>

struct SparseNode
{
  uint32_t parent;
  uint32_t children;
  Index global_index;
  uint4_t mask;
  // TODO Potential Neighbour caching
  // uint32_t begin;
  // uint32_t top;
  // uint32_t bottom;
  // uint32_t left;
  // uint32_t right;
};

struct TreeIndices
{

  uint8_t min_depth;
  uint8_t max_depth;
  SparseNode _indices[];

  uint32_t global_to_local(uint32_t global_index, int depth)
  {
    uint32_t index = 0;
    // iterate tree

    return index;
  };
};

template <typename T>
class QuadTree<T>
{
  TreeIndices& indices;
  uint64_t data_size;
  T _data[];
  T& operator[](Index I);
  // T& get_lowest_value(Index I); // get lowest depht value for global index
};

void broadcast_leafs(QuadTree<Index> index_cache, Range range)
{

  for (uint64_t i; i < index_cache.data_size; i++)
  {
    Index I = index_cache._data[i];
    // if Index in range:
    //       do something;
  }
};
void update_mipmap();

#endif // QUADTREE_H_
