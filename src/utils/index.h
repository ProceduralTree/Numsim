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

enum class Sign
{
  Plus,
  Minus
};

constexpr bool operator==(Offset lhs, Offset rhs)
{
  return lhs.x == rhs.x && lhs.y == rhs.y;
}
inline const Offset operator*(int n, const Offset& O) { return Offset { static_cast<int>(n * O.x), static_cast<int>(n * O.y) }; };
constexpr const Offset operator-(const Offset& O) { return { -O.x, -O.y }; };

struct Index
{
  uint16_t x;
  uint16_t y;
  uint16_t depth;

  inline Index operator+(const Index& other)
  {
    return { static_cast<uint16_t>(x + other.x), static_cast<uint16_t>(y + other.y), depth };
  };
  inline Index operator+(const Offset& other)
  {
    return { static_cast<uint16_t>(this->x + other.x), static_cast<uint16_t>(this->y + other.y), depth };
  };
  inline const Index operator+(const Offset& other) const { return { static_cast<uint16_t>(this->x + other.x), static_cast<uint16_t>(this->y + other.y) }; };
  inline Index operator-(const Offset& other)
  {
    assert(this->x >= other.x);
    assert(this->y >= other.y);
    return { static_cast<uint16_t>(this->x - other.x), static_cast<uint16_t>(this->y - other.y), depth };
  };
  inline const Index operator-(const Offset& other) const
  {
    assert(this->x >= other.x);
    assert(this->y >= other.y);
    return { static_cast<uint16_t>(this->x - other.x), static_cast<uint16_t>(this->y - other.y), depth };
  };
  inline bool operator<=(const Index& other) const
  {
    return x <= other.x && y <= other.y;
  }
  inline bool operator>=(const Index& other) const
  {
    return x >= other.x && y >= other.y;
  }
};

struct Range
{
  Index begin;
  Index end;
  inline Range operator+(const Offset& other)
  {
    return { begin + other, end + other };
  };
  inline size_t count() const
  {
    return (end.x - begin.x + 1) * (end.y - begin.y + 1);
  }
  inline Index size() const
  {
    return { static_cast<uint16_t>(end.x - begin.x + 1), static_cast<uint16_t>(end.y - begin.y + 1) };
  }
};

const inline constexpr Offset Ix = { 1, 0 };
const inline constexpr Offset Iy = { 0, 1 };
const inline constexpr Offset II { 1, 1 };

inline bool operator==(Index lhs, Index rhs)
{
  return lhs.x == rhs.x && lhs.y == rhs.y;
}
inline Range operator+(Range r, Offset o)
{
  return { r.begin + o, r.end + o };
}
inline Range operator-(Range r, Offset o)
{
  return { r.begin - o, r.end - o };
}

#endif // INDEX_H_
