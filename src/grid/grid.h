#ifndef GRID_H_
#define GRID_H_

#include <algorithm>
#include <cstdint>
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
    : top(Index { begin.x, end.y } + Iy, end + Iy)
    , bottom(begin - Iy, Index { end.x, begin.y } - Iy)
    , left(begin - Ix, Index { begin.x, end.y } - Ix)
    , right(Index { end.x, begin.y } + Ix, end + Ix)
    , all({ { top, Iy }, { bottom, -Iy }, { left, -Ix }, { right, Ix } })

  { };
  auto begin() const { return all.begin(); };
  auto end() const { return all.end(); };
};

class Grid2D
{
private:
  std::vector<double> _data;
  uint16_t size_x;
  uint16_t size_y;

public:
  Index begin;
  Index end;
  Offset len_x;
  Offset len_y;
  Range range;
  Boundaries boundary;

  Grid2D(uint16_t x, uint16_t y);
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
  double& operator[](Index index) { return (*this)[index.x, index.y]; };
  const double& operator[](Index index) const { return (*this)[index.x, index.y]; };
  double& operator[](uint16_t x, uint16_t y);
  const double& operator[](uint16_t x, uint16_t y) const;
  double& operator[](uint32_t z);
  const double& operator[](uint32_t z) const;

  const uint32_t elements() const { return this->_data.size(); }
  constexpr double max() { return *std::max_element(_data.begin(), _data.end()); }
  constexpr double min() { return *std::min_element(_data.begin(), _data.end()); }
};

void laplace(const Grid2D& in, Grid2D& out);
std::ostream& operator<<(std::ostream& os, const Grid2D& obj);

#endif // GRID_H_
