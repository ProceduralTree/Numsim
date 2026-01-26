#ifndef SPARSEGRID_H_
#define SPARSEGRID_H_
#include <cstddef>
#include <cstdint>
#include <grid/densetree.h>
#include <vector>

template <typename T>
struct SparseGrid2D
{
  const DenseTree::DenseTree& tree;
  std::vector<T> _data;

  constexpr double& operator[](size_t index)
  {
    return _data[index];
  };
  constexpr const double& operator[](size_t index) const
  {
    return _data[index];
  }

  constexpr double& operator[](Index I)
  {
    return _data[DenseTree::get_dense_index(tree, I)];
  };
  constexpr const double& operator[](Index I) const
  {
    return _data[DenseTree::get_dense_index(tree, I)];
  };
};

template <typename T>
void copy_entry(size_t index, uint16_t depth, const SparseGrid2D<T>& from, SparseGrid2D<T>& to)
{
  Zindex z = from._index_cache[index];
  size_t from_index = DenseTree::get_dense_index(to, z);
  to[index] = from[from_index];
}

template <typename T>
void copy(SparseGrid2D<T> from, SparseGrid2D<T> to)
{
  DenseTree::broadcast_breath_first(copy_entry, to.tree, to.tree.max, from, to);
}

#endif // SPARSEGRID_H_
