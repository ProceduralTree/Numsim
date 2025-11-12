#ifndef INDEX_H_
#define INDEX_H_

#include <cassert>
#include <cstdint>
#include <functional>

struct Offset
{
  int x;
  int y;
};

constexpr bool operator==(Offset lhs, Offset rhs)
{
  return lhs.x == rhs.x && lhs.y == rhs.y;
}
constexpr const Offset operator*(int n, const Offset& O) { return Offset { static_cast<int>(n * O.x), static_cast<int>(n * O.y) }; };
constexpr const Offset operator-(const Offset& O) { return -1 * O; };

struct Index
{
  uint16_t x;
  uint16_t y;

  constexpr Index operator+(const Offset& other)
  {
    return { static_cast<uint16_t>(this->x + other.x), static_cast<uint16_t>(this->y + other.y) };
  };
  constexpr const Index operator+(const Offset& other) const { return { static_cast<uint16_t>(this->x + other.x), static_cast<uint16_t>(this->y + other.y) }; };
  constexpr Index operator-(const Offset& other)
  {
    assert(this->x >= other.x);
    assert(this->y >= other.y);
    return { static_cast<uint16_t>(this->x - other.x), static_cast<uint16_t>(this->y - other.y) };
  };
  constexpr const Index operator-(const Offset& other) const
  {
    assert(this->x >= other.x);
    assert(this->y >= other.y);
    return { static_cast<uint16_t>(this->x - other.x), static_cast<uint16_t>(this->y - other.y) };
  };
};

struct Range
{
  Index begin;
  Index end;
  constexpr Range operator+(const Offset& other)
  {
    return { begin + other, end + other };
  };
};

constexpr Offset Ix = { 1, 0 };
constexpr Offset Iy = { 0, 1 };

constexpr bool operator==(Index lhs, Index rhs)
{
  return lhs.x == rhs.x && lhs.y == rhs.y;
}

#endif // INDEX_H_
