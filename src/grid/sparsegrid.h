#ifndef SPARSEGRID_H_
#define SPARSEGRID_H_
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <grid/densetree.h>
#include <iterator>
#include <vector>

template <typename T>
struct SparseGrid2D
{
  const DenseTree::DenseTree& tree;
  std::vector<T> _data;

  constexpr T& operator[](size_t index)
  {
    return _data.at(index);
  };
  constexpr const T& operator[](size_t index) const
  {
    return _data.at(index);
  }

  constexpr T& operator[](Index I)
  {
    return _data.at(DenseTree::get_dense_index(tree, I));
  };
  constexpr const T& operator[](Index I) const
  {
    return _data.at(DenseTree::get_dense_index(tree, I));
  }
  SparseGrid2D(const SparseGrid2D&) = delete;
  SparseGrid2D& operator=(const SparseGrid2D&) = delete;

  SparseGrid2D(SparseGrid2D&&) = default;
  SparseGrid2D& operator=(SparseGrid2D&&) = default;

  SparseGrid2D(const DenseTree::DenseTree& tree)
    : tree(tree)
    , _data(tree._sizes[tree.maxDepth + 1], 0) { };

  SparseGrid2D(const DenseTree::DenseTree& tree, std::vector<T> data)
    : tree(tree)
    , _data(tree._sizes[tree.maxDepth + 1], 0)
  {
    assert(data.size() == _data.size());
    std::copy(data.begin(), data.end(), _data);
  };
  constexpr friend void swap(SparseGrid2D& a, SparseGrid2D& b) noexcept
  {
    assert(&a.tree == &b.tree && "SparseGrid2D swap: tree references must match");
    using std::swap;
    swap(a._data, b._data);
  }
};

template <typename T>
void copy_entry(size_t index, uint16_t depth, const SparseGrid2D<T>& from, SparseGrid2D<T>& to)
{
  Zindex z = from._index_cache[index];
  size_t from_index = DenseTree::get_dense_index(to, z);
  to[index] = from[from_index];
}

template <typename Operator, typename T, typename... Args>
void mipmap(Operator&& O, SparseGrid2D<T>& grid, Args&&... args)
{

  for (uint16_t depth = grid.tree.maxDepth; depth > 0; depth--)
  {
    // std::cerr << "Depth:" << depth << std::endl;
    for (size_t local_index = grid.tree._sizes.at(depth - 1); local_index < grid.tree._sizes.at(depth); local_index++)
    {
      if (grid.tree._depths[local_index] > 0)
      {
        size_t data_index = grid.tree._indices[local_index];
        grid[local_index] = std::forward<Operator>(O)({ grid[data_index], grid[data_index + 1], grid[data_index + 2], grid[data_index + 3] }, std::forward<Args>(args)...);
      }
    }
  }
}

template <typename T>
T _sum(std::array<T, 4> data)
{
  return data[0] + data[1] + data[2] + data[3];
};
template <typename T>
T _mean(std::array<T, 4> data)
{
  return 0.25 * (data[0] + data[1] + data[2] + data[3]);
};
template <typename T>
T _max(std::array<T, 4> data)
{
  return std::max(data[0], data[1], data[2], data[3]);
};
template <typename T>
T _or(std::array<T, 4> data)
{
  return (data[0] | data[1] | data[2] | data[3]);
};

// template <typename T>
// void copy(SparseGrid2D<T> from, SparseGrid2D<T> to)
//{
//   DenseTree::broadcast_breath_first(copy_entry, to.tree, to.tree.max, from, to);
// }

#endif // SPARSEGRID_H_
