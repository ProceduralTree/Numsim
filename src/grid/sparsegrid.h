#ifndef SPARSEGRID_H_
#define SPARSEGRID_H_
#include <grid/densetree.h>

class SparseGrid2D
{
  const DenseTree::DenseTree& _indices;
  double _data[];

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
    return _data[DenseTree::get_dense_index(_indices, I)];
  };
  constexpr const double& operator[](Index I) const
  {
    return _data[DenseTree::get_dense_index(_indices, I)];
  };
};

#endif // SPARSEGRID_H_
