#ifndef QUADTREE_H_
#define QUADTREE_H_
#include "utils/index.h"
#include "utils/stb_image.h"
#include <bits/floatn.h>
#include <cstdint>
#include <filesystem>
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
    size_t oldSize = _numBytes;
    _numBits = bitSize;
    if (bitSize < 8)
      _numBytes = 1;
    else
      _numBytes = 1 + (bitSize - 1) / 8;
    _data = (uint8_t*)realloc(_data, _numBytes);
    std::fill(_data + oldSize, _data + _numBytes, 0);
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
  size_t _numBits = 0;
  size_t _numBytes = 0;
};

inline bool isPowerOf2(size_t x)
{
  return __builtin_popcount(x) == 1;
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
  QuadTree<T>(std::filesystem::path imagePath)
  {
    generateGridFromImage(imagePath);
  }
  QuadTree<T>(size_t sizeX, size_t sizeY)
  {
    createWithSize(sizeX, sizeY);
  }
  ~QuadTree<T>()
  {
    free(_data);
  }

  const T& getLowestValue(Index I) const
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    size_t d = depth - 1;
    for (; d > 1 && !tree[(zorder >> (2 * d)) + depthOffset[d]]; d--)
    {
    }
    return _data[(zorder >> (2 * d)) + depthOffset[d]];
  }
  const size_t getLowestDepth(Index I) const
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    size_t d = depth - 1;
    for (; d > 1 && !tree[(zorder >> (2 * d)) + depthOffset[d]]; d--)
    {
    }
    return d;
  }
  T& operator[](Index I)
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    return _data[(zorder >> (2 * I.depth)) + depthOffset[I.depth]];
  }
  const T& operator[](Index I) const
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    return _data[(zorder >> (2 * I.depth)) + depthOffset[I.depth]];
  }
  void SetTreeNode(Index I, bool value = true)
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    tree.set((zorder >> (2 * I.depth)) + depthOffset[I.depth], value);
  }
  const bool GetTreeNode(Index I) const
  {
    size_t zorder = IndexToZOrder(I.x, I.y);
    return tree.get((zorder >> (2 * I.depth)) + depthOffset[I.depth]);
  }

private:
  dynamic_bitset tree;
  size_t depth;
  T* _data;
  std::vector<size_t> depthOffset;
  size_t maxSize;

  size_t clamp(size_t min, size_t max, size_t x)
  {
    return min > x ? min : max < x ? max
                                   : x;
  }
  bool generateGridFromImage(std::filesystem::path imagePath, int rgbOffset = 0)
  {

    int width, height, channels;
    unsigned char* data = stbi_load(imagePath.c_str(), &width, &height, &channels, 4);
    size_t leadingZeros = std::min(__builtin_clz(width), __builtin_clz(height)); // ggrks
    size_t size = ((((size_t)-1) >> 1) + 1) >> (leadingZeros - 1);
    createWithSize(size, size);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < width; y++)
      {
        int8_t down = (int8_t)data[(width * clamp(0, height - 1, y - 1) + x) * channels]; // this is boundary check down
        int8_t left = (int8_t)data[(width * y + clamp(0, width - 1, x - 1)) * channels]; // this is boundary check left
        int8_t d = (int8_t)data[(width * y + x) * channels + rgbOffset];
        if (rgbOffset == 0) // reading p
        {
          if (d != 0)
            SetTreeNode({ x, y, depth });
        } else if (rgbOffset == 1) // reading u
        {
          if (d != 0)
          {
            if (left != 0)
            {
              SetTreeNode({ x + 1, y, depth });
              operator[]({ x + 1, y, depth }) = data;

            } else
            {
              SetTreeNode({ x, y, depth });
              operator[]({ x, y, depth }) = data;
            }
          }

        } else if (rgbOffset == 2) // reading v
        {
          if (d != 0)
          {
            if (down != 0)
            {
              SetTreeNode({ x, y + 1, depth });
              operator[]({ x, y + 1, depth }) = data;

            } else
            {
              SetTreeNode({ x, y, depth });
              operator[]({ x, y, depth }) = data;
            }
          }
        }
      }
    }
    floodTreeToRoot();
  }
  void createWithSize(size_t sizeX, size_t sizeY)
  {
    assert(isPowerOf2(sizeX) && sizeX == sizeY);
    maxSize = sizeX;
    depth = __builtin_ctz(sizeX);
    size_t allocSize = sizeX * sizeY;
    size_t tempSizeX = sizeX / 2;
    size_t tempSizeY = sizeY / 2;
    for (size_t i = 1; i < depth - 1; ++i) // from 1.. d-1 because max depth is already in allocSize and root node is 4 not 1
    {
      allocSize += tempSizeX * tempSizeY;
      tempSizeX /= 2;
      tempSizeY /= 2;
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
  void floodTreeToRoot()
  {
    size_t offset = 1;
    for (size_t d = depth - 1; d >= 0; d--)
    {
      offset *= 2;
      for (size_t x = 0; x < maxSize; x += offset)
      {
        for (size_t y = 0; y < maxSize; y += offset)
        {
          Index I = { (uint16_t)x, (uint16_t)y, (uint8_t)(d + 1) };
          size_t zorder = IndexToZOrder(I.x, I.y);
          size_t index = (zorder >> (2 * I.depth)) + depthOffset[I.depth];
          if (tree[index] | tree[index + 1] | tree[index + 2] | tree[index + 3])
          {
            SetTreeNode({ x, y, d });
          }
        }
      }
    }
  }
};

#endif // QUADTREE_H_
