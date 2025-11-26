#ifndef GRID_H_
#define GRID_H_

#include <algorithm>
#include <cstdint>
#include <execution>
#include <iostream>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "indexing.h"
#include "utils/index.h"

#define CARTESIAN

struct Boundaries
{
  const Range top;
  const Range bottom;
  const Range left;
  const Range right;
  const std::array<std::tuple<Range, Offset>, 4> all;
  Boundaries(Index begin, Index end)
    : top(Index { begin.x, end.y } + Iy - Ix, end + Iy + Ix)
    , bottom(begin - Iy - Ix, Index { end.x, begin.y } + Ix - Iy)
    , left(begin - Ix - Iy, Index { begin.x, end.y } - Ix + Iy)
    , right(Index { end.x, begin.y } + Ix - Iy, end + Ix + Iy)
    , all({ { top, Iy }, { bottom, -Iy }, { left, -Ix }, { right, Ix } })

  { };
  auto begin() const { return all.begin(); };
  auto end() const { return all.end(); };
};

class Grid2D
{

public:
  uint16_t size_x;
  uint16_t size_y;
  Index begin;
  Index end;
  Range range;
  Boundaries boundary;

  Grid2D(Index beg, Index end);

  Grid2D(const Grid2D&) = delete;
  Grid2D& operator=(const Grid2D&) = delete;

  Grid2D(Grid2D&& other) noexcept
    : boundary(other.boundary)
  {
    std::swap(size_x, other.size_x);
    std::swap(size_y, other.size_y);
    std::swap(_data, other._data);
  }

  Grid2D& operator=(Grid2D&& other) noexcept
  {
    std::swap(size_x, other.size_x);
    std::swap(size_y, other.size_y);
    std::swap(_data, other._data);
    return *this;
  }

  // double& operator[](Index index) { return (*this)[index.x + begin.x, index.y + begin.y]; };
  // const double& operator[](Index index) const { return (*this)[index.x + begin.x, index.y + begin.y]; };
  inline double& operator[](Index I)
  {
#ifdef CARTESIAN
    uint32_t index = I.x + size_x * I.y;
#else
    uint32_t index = z_order(x, y);
#endif
#ifdef DEBUG
    return this->_data.at(index);
#else
    return this->_data.data()[index];
#endif
  };

  inline const double& operator[](Index I) const
  {
#ifdef CARTESIAN
    uint32_t index = I.x + size_x * I.y;
#else
    uint32_t index = z_order(x, y);
#endif
#ifdef DEBUG
    return this->_data.at(index);
#else
    return this->_data.data()[index];
#endif
  };

  double& operator[](uint32_t z);
  const double& operator[](uint32_t z) const;

  const uint32_t elements() const { return this->_data.size(); }
  inline double max()
  {
    return *std::max_element(std::execution::par_unseq, _data.begin(), _data.end());
    // double result = 0;

    // #pragma omp parallel for collapse(2) reduction(max : result)
    // for (uint16_t j = begin.y; j <= end.y; j++)
    //{
    // for (uint16_t i = begin.x; i <= end.x; i++)
    //{
    // Index I = { i, j };
    // result = std::max(result, (*this)[I]);
    //}
    //}
    // return result;
  };
  inline double min()
  {
    return *std::min_element(std::execution ::par_unseq, _data.begin(), _data.end());
    // double result = 0;

    // #pragma omp parallel for collapse(2) reduction(min : result)
    // for (uint16_t j = begin.y; j <= end.y; j++)
    //{
    // for (uint16_t i = begin.x; i <= end.x; i++)
    //{
    // Index I = { i, j };
    // result = std::max(result, (*this)[I]);
    // }
    // }
    // return result;
  };

private:
  std::vector<double> _data;
};
std::ostream& operator<<(std::ostream& os, const Grid2D& obj);

#endif // GRID_H_
