#ifndef INDEXING_H_
#define INDEXING_H_

#include <cstdint>
#include <tuple>

constexpr uint32_t masky = 0xAAAA;
constexpr uint32_t maskx = 0x5555;

enum class Indexing
{
  Cartesian,
  ZOrder,
};

struct Indices
{
  uint32_t top;
  uint32_t bottom;
  uint32_t left;
  uint32_t right;
};

// Spread 16 bits apart by inserting zeros between bits (bitwise magic)
constexpr uint32_t part1by1(uint16_t n)
{
  uint32_t x = n;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
}

// Morton encode: interleave bits of x and y
constexpr uint32_t z_order(uint16_t x, uint16_t y)
{
  return (part1by1(y) << 1) | part1by1(x);
}

// Compact every other bit into a single 16-bit number
constexpr uint16_t compact1by1(uint32_t n)
{
  n &= 0x55555555;
  n = (n | (n >> 1)) & 0x33333333;
  n = (n | (n >> 2)) & 0x0F0F0F0F;
  n = (n | (n >> 4)) & 0x00FF00FF;
  n = (n | (n >> 8)) & 0x0000FFFF;
  return static_cast<uint16_t>(n);
};

// Morton decode
constexpr std::pair<uint16_t, uint16_t> decode_z_order(uint32_t z)
{
  uint16_t x = compact1by1(z);
  uint16_t y = compact1by1(z >> 1);
  return { x, y };
};

constexpr Indices indices(uint32_t z)
{
  uint32_t x_bits_masked = z & maskx;
  uint32_t y_bits_masked = z & masky;
  uint32_t top = (x_bits_masked - 1) & maskx;
  uint32_t bottom = (y_bits_masked + 1) & maskx;
  uint32_t left = (y_bits_masked - 1) & masky;
  uint32_t right = (x_bits_masked + 1) & masky;
  top |= y_bits_masked;
  bottom |= y_bits_masked;
  left |= x_bits_masked;
  right |= x_bits_masked;
  return { top, bottom, left, right };
};

uint32_t interleave(uint16_t x, uint16_t y);
std::tuple<uint32_t, uint32_t> detangle(uint32_t z);
#endif // INDEXING_H_
