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

#define CARTESIAN

class Grid2D
{
private:
  std::vector<double> _data;
  uint16_t size_x;
  uint16_t size_y;

public:
  Grid2D(uint16_t x, uint16_t y);

  Grid2D(const Grid2D&) = delete;
  Grid2D& operator=(const Grid2D&) = delete;

  Grid2D(Grid2D&& other) noexcept
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

  double& operator[](uint16_t x, uint16_t y);
  const double& operator[](uint16_t x, uint16_t y) const;
  double& operator[](uint32_t z);
  const double& operator[](uint32_t z) const;

  const uint32_t elements() const { return this->_data.size(); }
  constexpr double max() { return *std::max_element(_data.begin(), _data.end()); }
};

void laplace(const Grid2D& in, Grid2D& out);
std::ostream& operator<<(std::ostream& os, const Grid2D& obj);

#endif // GRID_H_
