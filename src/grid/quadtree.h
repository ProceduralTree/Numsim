#ifndef QUADTREE_H_
#define QUADTREE_H_
#include <cstdint>

struct SparseNode
{
  uint32_t parent;
  uint32_t children;
  uint4_t mask;
  // TODO Potential Neighbour caching
  // uint32_t begin;
  // uint32_t top;
  // uint32_t bottom;
  // uint32_t left;
  // uint32_t right;
};

struct treeIndices
{
  uint8_t levels;
  SparseNode _indices[];

  uint32_t global_to_local(uint32_t global_index)
  {
  }
};

template <typename T>
class QuadTree<T>
{
  T _data[];
};

#endif // QUADTREE_H_
