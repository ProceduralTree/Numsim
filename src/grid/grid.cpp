#include "grid.h"
#include "utils/index.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>

Grid2D::Grid2D(Index beg, Index end)
  : _data((end.x + 2) * (end.y + 2), 0.)
  , size_x(end.x - beg.x + 3)
  , size_y(end.y - beg.y + 3)
  , begin(beg)
  , end(end)
  , range(beg, end)
  , boundary(beg, end)
{
  // uint32_t size = std::bit_width(x) + std::bit_width(y);
  // this->_data.resize(1 << size, 0.);
  // this->_data.resize(x * y, init);
}

double& Grid2D::operator[](uint32_t index) { return this->_data[index]; }

const double& Grid2D::operator[](uint32_t index) const
{
#ifdef NDEBUG
  return this->_data[index];
#else
  return this->_data.at(index);
#endif
}

std::ostream& operator<<(std::ostream& os, const Grid2D& obj)
{
  os << std::scientific << std::setprecision(3) << std::endl;
  os << (obj.end.x - obj.begin.x) << "x" << (obj.end.y - obj.begin.y) << " Grid2D" << std::endl;

  const int width = 5;
  const int len = 10;

  for (uint16_t j = obj.begin.y - 1; j < std::min(obj.begin.y + width, static_cast<int>(obj.end.y)); j++)
  {
    for (uint16_t i = obj.begin.x - 1; i < std::min(obj.begin.x + width, static_cast<int>(obj.end.x)); i++)
    {
      os << std::setw(len) << obj[Index { i, j }] << "";
    }
    if (obj.end.x - obj.begin.x > 2 * width + 2)
    {
      os << std::setw(3) << "  …";
    }
    for (uint16_t i = std::max(obj.end.x - width, obj.begin.x + width); i < obj.end.x + 2; i++)
    {
      os << std::setw(len) << obj[Index { i, j }] << "";
    }
    os << std::endl;
  }
  if (obj.end.y - obj.begin.y > 2 * width + 2)
  {
    for (uint16_t i = obj.begin.x - 1; i < std::min(obj.begin.x + width, static_cast<int>(obj.end.x)); i++)
    {
      os << std::setw(len) << "  ⋮" << "";
    }
    if (obj.end.x - obj.begin.x > 2 * width + 2)
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
      os << std::setw(len) << obj[{ i, j }] << "";
    }
    if (obj.end.x - obj.begin.x > 2 * width + 2)
    {
      os << std::setw(3) << "  …";
    }
    for (uint16_t i = std::max(obj.end.x - width, obj.begin.x + width); i < obj.end.x + 2; i++)
    {
      os << std::setw(len) << obj[{ i, j }] << "";
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
