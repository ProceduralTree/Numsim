#ifndef GRID_H_
#define GRID_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <execution>
#include <iostream>
#include <mpi.h>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "indexing.h"
#include "utils/Logger.h"
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
  inline std::array<std::tuple<Range, Offset>, 4> u_ghosts()
  {
    return { std::tuple<Range, Offset> { top, Iy }, { bottom, -Iy }, { left, -2 * Ix }, { right, 2 * Ix } };
  };
  inline std::array<std::tuple<Range, Offset>, 4> v_ghosts()
  {
    return { std::tuple<Range, Offset> { top, 2 * Iy }, { bottom, -2 * Iy }, { left, -Ix }, { right, Ix } };
  };
  inline std::array<std::tuple<Range, Offset>, 4> unique()
  {
    return { std::tuple<Range, Offset> { Range { top.begin + Ix, top.end - Ix }, Iy }, { { bottom.begin + Ix, bottom.end - Ix }, Iy }, { left, -Ix }, { right, Ix } };
  };
};

class Grid2D
{

public:
  uint16_t size_x;
  uint16_t size_y;
  Index begin;
  Index end;
  Range range;
  Range globalRange;
  Boundaries boundary;

  Grid2D(Index beg, Index end);
  Grid2D(Index beg, Index end, Range globalRange);

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
    assert(I <= range.end + II && "invalid grid range");
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
    if (!(I <= range.end + II))
    {
      DebugF("invalid grid range size:{},{}, index {},{}", size_x, size_y, I.x, I.y);
    }
    assert(I <= range.end + II && "invalid grid range");
    return this->_data.at(index);
#else
    return this->_data.data()[index];
#endif
  };

  double& operator[](uint32_t z);
  const double& operator[](uint32_t z) const;

  void get(double* buffer, Range r) const
  {
    assert(r.begin.x >= begin.x - 1);
    assert(r.begin.y >= begin.y - 1);
    assert(r.end.x <= end.x + 1);
    assert(r.end.y <= end.y + 1);

    uint32_t index = 0;
    for (uint16_t j = r.begin.y; j <= r.end.y; j++)
    {
      for (uint16_t i = r.begin.x; i <= r.end.x; i++, index++)
      {
        buffer[index] = (*this)[{ i, j }];
      }
    }
  };
  inline void set(double* buffer, Range r)
  {
    assert(r.begin.x >= begin.x - 1);
    assert(r.begin.y >= begin.y - 1);
    assert(r.end.x <= end.x + 1);
    // DebugF("end {{x={},y={}}}  , r {{x={},y={}}}", end.x, end.y, r.end.x, r.end.y);
    assert(r.end.y <= end.y + 1);

    uint32_t index = 0;
    for (uint16_t j = r.begin.y; j <= r.end.y; j++)
    {
      for (uint16_t i = r.begin.x; i <= r.end.x; i++, index++)
      {
        (*this)[{ i, j }] = buffer[index];
      }
    }
  };

  const uint32_t elements() const { return this->_data.size(); }
  inline double max()
  {
    double local_max = *std::max_element(_data.begin(), _data.end());
    double global_max = 0.;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_max;
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
    double local_min = *std::min_element(_data.begin(), _data.end());
    double global_min = 0.;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global_min;
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
