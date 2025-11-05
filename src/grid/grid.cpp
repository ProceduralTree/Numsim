#include "grid.h"
#include "indexing.h"
#include "utils/broadcast.h"
#include "utils/index.h"
#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>

Grid2D::Grid2D(uint16_t x, uint16_t y)
  : _data(x * y, 0.)
  , size_x(x)
  , size_y(y)
  , begin(0, 0)
  , end(size_x, size_y)
  , len_x(size_x, 0)
  , len_y(size_y, 0)
  , range(begin, end) {
    // uint32_t size = std::bit_width(x) + std::bit_width(y);
    // this->_data.resize(1 << size, 0.);
    // this->_data.resize(x * y, init);
  };
Grid2D::Grid2D(Index beg, Index end)
  : _data((end.x + 2) * (end.y + 2), 0.)
  , size_x(end.x - beg.x + 3)
  , size_y(end.y - beg.y + 3)
  , begin(beg)
  , end(end)
  , len_x(size_x - 3, 0)
  , len_y(0, size_y - 3)
  , range(beg, end) {
    // uint32_t size = std::bit_width(x) + std::bit_width(y);
    // this->_data.resize(1 << size, 0.);
    // this->_data.resize(x * y, init);
  };

double& Grid2D::operator[](uint16_t x, uint16_t y)
{
#ifdef CARTESIAN
  uint32_t index = x + size_x * y;
#else
  std::cout << "!ZORDER" << std::endl;
  uint32_t index = z_order(x, y);
#endif
  return this->_data.at(index);
}

const double& Grid2D::operator[](uint16_t x, uint16_t y) const
{
#ifdef CARTESIAN
  uint32_t index = x + size_x * y;
#else
  uint32_t index = z_order(x, y);
#endif
  return this->_data.at(index);
}

double& Grid2D::operator[](uint32_t index) { return this->_data[index]; };

const double& Grid2D::operator[](uint32_t index) const
{
  return this->_data.at(index);
};

std::ostream& operator<<(std::ostream& os, const Grid2D& obj)
{
  // std::cout << "\n";
  // for (int j = obj.begin.y - 1; j < obj.end.y + 2; j++)
  //{
  //   for (int i = obj.begin.x - 1; i < obj.end.x + 2; i++)
  //   {
  //     std::cout << obj[i, j] << "\t";
  //   }
  //   std::cout << "\n";
  // }
  //
  os << std::scientific << std::setprecision(1) << std::endl;
  os << obj.size_x << "x" << obj.size_y << " Grid2D" << std::endl;

  const int width = 5;
  const int len = 10;

  for (uint16_t j = obj.begin.y - 1; j < std::min(obj.begin.y + width, static_cast<int>(obj.end.y)); j++)
  {
    for (uint16_t i = obj.begin.x - 1; i < std::min(obj.begin.x + width, static_cast<int>(obj.end.x)); i++)
    {
      os << std::setw(len) << obj[i, j] << "";
    }
    if (obj.len_x.x > 2 * width + 2)
    {
      os << std::setw(3) << "  …";
    }
    for (uint16_t i = std::max(obj.end.x - width, obj.begin.x + width); i < obj.end.x + 2; i++)
    {
      os << std::setw(len) << obj[i, j] << "";
    }
    os << std::endl;
  }
  if (obj.len_y.y > 2 * width + 2)
  {
    for (uint16_t i = obj.begin.x - 1; i < std::min(obj.begin.x + width, static_cast<int>(obj.end.x)); i++)
    {
      os << std::setw(len) << "  ⋮" << "";
    }
    if (obj.len_x.x > 2 * width + 2)
    {
      os << std::setw(3) << "  ⋱";
    }
    for (uint16_t i = std::max(obj.end.x - width, obj.begin.x + width); i < obj.end.x + 2; i++)
    {
      os << std::setw(len) << "  ⋮" << "";
    }
    os << std::endl;
  }
  for (uint16_t j = std::max(obj.end.y - width, obj.begin.y + width); j < obj.end.y + 2; j++)
  {
    for (uint16_t i = obj.begin.x - 1; i < std::min(obj.begin.x + width, static_cast<int>(obj.end.x)); i++)
    {
      os << std::setw(len) << obj[i, j] << "";
    }
    if (obj.len_x.x > 2 * width + 2)
    {
      os << std::setw(3) << "  …";
    }
    for (uint16_t i = std::max(obj.end.x - width, obj.begin.x + width); i < obj.end.x + 2; i++)
    {
      os << std::setw(len) << obj[i, j] << "";
    }
    os << std::endl;
  }

  return os;
}
// bool boundary(uint32_t zindex, uint16_t sx, uint16_t sy)
//{
//   auto [x, y] = decode_z_order(zindex);
//   return ((x < sx) && (x > 0)) || ((y < sy) && (y > 0));
// }
