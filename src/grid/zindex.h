#ifndef ZIDEX_H_
#define ZIDEX_H_

#include <utils/index.h>

#include <cstddef>
#include <cstdint>

constexpr size_t part1by1(const size_t n)
{
  size_t x = n;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
};

constexpr size_t IndexToZOrder(const size_t x, const size_t y)
{
  return (part1by1(y) << 1) | part1by1(x);
};

inline uint16_t compact1by1(uint32_t n)
{
  n &= 0x55555555;
  n = (n ^ (n >> 1)) & 0x33333333;
  n = (n ^ (n >> 2)) & 0x0F0F0F0F;
  n = (n ^ (n >> 4)) & 0x00FF00FF;
  n = (n ^ (n >> 8)) & 0x0000FFFF;
  return static_cast<uint16_t>(n);
};

constexpr Index ZorderToIndex(uint32_t index)
{
  uint16_t x = compact1by1(index >> 0); // even bits
  uint16_t y = compact1by1(index >> 1); // odd bits
  return Index { x, y, 0 };
};

struct Zindex
{
  size_t index;
  uint16_t depth;
  Zindex(Index I)
    : index(IndexToZOrder(I.x, I.y))
    , depth(I.depth) { };
  Zindex(size_t z, uint16_t d)
    : index(z)
    , depth(d) { };
};

#endif // ZIDEX_H_
