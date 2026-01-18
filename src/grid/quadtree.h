#ifndef QUADTREE_H_
#define QUADTREE_H_
#include "utils/index.h"
#include <cstdint>
#include <sstream>

struct dynamic_bitset
{
  dynamic_bitset() { }
  dynamic_bitset(const size_t bitSize)
  {
    resize(bitSize);
  }
  ~dynamic_bitset()
  {
    free(_data);
  }

  void set(const size_t index, const bool value = true)
  {
    size_t byteIndex = index / 8;
    size_t offset = index % 8;

    uint8_t bitfield = 1 << offset;
    if (value)
      _data[byteIndex] |= bitfield;
    else
      _data[byteIndex] &= ~bitfield;
  }
  void setAll(const bool value)
  {
    uint8_t fillValue = 0xFF;
    if (!value)
      fillValue = 0;
    std::fill(_data, _data + _numBytes, fillValue);
  }
  void reset()
  {
    setAll(false);
  }
  const bool get(const size_t index) const
  {
    assert(index < _numBits);
    return operator[](index);
  }
  void flip(const size_t index)
  {
    set(index, !operator[](index));
  }

  const bool operator[](size_t index) const
  {
    size_t byteIndex = index / 8;
    size_t offset = index % 8;
    return (_data[byteIndex] >> offset) & 0x1;
  }

  void resize(size_t bitSize)
  {
    _numBits = bitSize;
    if (bitSize < 8)
      _numBytes = 1;
    else
      _numBytes = 1 + (bitSize - 1) / 8;
    _data = (uint8_t*)realloc(_data, _numBytes);
  }
  const uint8_t* data() const
  {
    return _data;
  }
  const size_t size() const
  {
    return _numBits;
  }
  const size_t dataSize() const
  {
    return _numBytes;
  }

  std::string toString() const
  {
    std::stringstream ss;
    for (size_t i = 0; i < _numBits; ++i)
    {
      ss << (operator[](i) ? "1" : "0");
    }
    return ss.str();
  }

private:
  uint8_t* _data = nullptr;
  size_t _numBits;
  size_t _numBytes;
};

inline bool isPowerOf4(size_t x)
{
  return x != 0 && ((x & (x - 1)) == 0) && !(x & 0xAAAAAAAA);
}
inline size_t part1by1(const size_t n)
{
  size_t x = n;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
}

inline size_t IndexToZOrder(const size_t x, const size_t y)
{
  return (part1by1(y) << 1) | part1by1(x);
}

template <typename T>
struct QuadTree
{
  QuadTree<T>(size_t sizeX, size_t sizeY)
  {
    assert(isPowerOf4(sizeX) && sizeX == sizeY);
    depth = __builtin_ctz(sizeX) / 2;
    size_t allocSize = sizeX * sizeY;
    size_t tempSizeX = sizeX / 4;
    size_t tempSizeY = sizeY / 4;
    for (size_t i = 1; i < depth; ++i)
    {
      allocSize += tempSizeX * tempSizeY;
      tempSizeX = sizeX / 4;
      tempSizeY = sizeY / 4;
    }

    _data = malloc(sizeof(T) * allocSize);
    tree.resize(allocSize - sizeX * sizeY);
    depthOffset.resize(depth);
    size_t offset = 0;
    size_t nodesCount = 4;
    for (size_t i = 0; i < depth; i++)
    {
      depthOffset[i] = offset;
      offset += nodesCount;
      nodesCount *= 4;
    }
  }
  ~QuadTree<T>()
  {
    free(_data);
  }

  size_t calcTreeIndex(Index I)
  {
    size_t DepthOffset = I.depth * 4;
  }
  T& operator[](Index I)
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    size_t d = depth - 1;
    for (size_t d = depth - 1; d > 1 && !tree[(zorder >> 2 * d) + depthOffset[d]]; d--)
    {
    }
    return _data[(zorder >> 2 * d) + depthOffset[d]];
  }
  const T& operator[](Index I) const
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    size_t d = depth - 1;
    for (size_t d = depth - 1; d > 1 && !tree[(zorder >> 2 * d) + depthOffset[d]]; d--)
    {
    }
    return _data[(zorder >> 2 * d) + depthOffset[d]];
  }

private:
  dynamic_bitset tree;
  size_t depth;
  T* _data;
  std::vector<size_t> depthOffset;
};

#endif // QUADTREE_H_
