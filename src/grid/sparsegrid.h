#ifndef SPARSEGRID_H_
#define SPARSEGRID_H_
#include <grid/densetree.h>

class SparseGrid2D
{
  const DenseTree::DenseTree& _indices;
  double _data[];

  constexpr double& operator[](Index I)
  {
    DenseTree::get_dense_index(_indices, idx);
  };
  constexpr const double& operator[](Index I) const;
};

#endif // SPARSEGRID_H_
