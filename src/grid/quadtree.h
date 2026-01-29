#ifndef QUADTREE_H_
#define QUADTREE_H_
#include "utils/Logger.h"
#include "utils/index.h"
#include "utils/stb_image.h"
#include <bits/floatn.h>
#include <cstdint>
#include <filesystem>
#include <iostream>
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
    _data = nullptr;
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
inline size_t part1by1_l(const size_t n)
{
  size_t x = n;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
}

inline size_t IndexToZOrder_l(const size_t x, const size_t y)
{
  return (part1by1_l(y) << 1) | part1by1_l(x);
}

template <typename T>
struct QuadTree
{
  QuadTree<T>(const std::filesystem::path& imagePath)
  {
    generateGridFromImage(imagePath);
  }
  QuadTree<T>(size_t sizeX, size_t sizeY)
  {
    createWithSize(sizeX, sizeY);
  }
  ~QuadTree<T>()
  {
    free(_dataU);
    _dataU = nullptr;
    free(_dataV);
    _dataV = nullptr;
  }
  QuadTree<T>(const QuadTree<T>&) = delete;
  QuadTree<T>& operator=(const QuadTree<T>&) = delete;
  QuadTree<T>(QuadTree<T>&&) = delete;
  QuadTree<T>& operator=(QuadTree<T>&&) = delete;

  size_t getZorderIndex(Index I) const
  {
    size_t zorder = IndexToZOrder_l(I.x, I.y);
    return (zorder >> (2 * (depth - I.depth))) + depthOffset[I.depth - 1];
  }
  const T& getULowestValue(Index I) const
  {
    I.depth = 0;
    while (I.depth < depth && uTree.get(getZorderIndex(I)))
    {
      I.depth++;
    }
    return _dataU[getZorderIndex(I)];
  }
  const size_t getPLowestDepth(Index I) const
  {
    I.depth = 0;
    while (I.depth < depth && pTree.get(getZorderIndex(I)))
    {
      I.depth++;
    }
    return I.depth;
  }

  void SetPTreeNode(Index I, bool value = true)
  {
    size_t index = getZorderIndex(I);
    pTree.set(index, value);
  }
  const bool GetPTreeNode(Index I) const
  {
    size_t index = getZorderIndex(I);
    return pTree.get(index);
  }
  void SetUTreeNode(Index I, bool value = true)
  {
    size_t index = getZorderIndex(I);
    uTree.set(index, value);
  }
  const bool GetUTreeNode(Index I) const
  {
    size_t index = getZorderIndex(I);
    return uTree.get(index);
  }
  void SetVTreeNode(Index I, bool value = true)
  {
    size_t index = getZorderIndex(I);
    vTree.set(index, value);
  }
  const bool GetVTreeNode(Index I) const
  {
    size_t index = getZorderIndex(I);
    return vTree.get(index);
  }
  void SetPData(Index I, bool value = true)
  {
    size_t index = getZorderIndex(I);
    return _dataP.set(index, value);
  }
  void SetUData(Index I, T data)
  {
    size_t index = getZorderIndex(I);
    _dataU[index] = data;
  }
  void SetVData(Index I, T data)
  {
    size_t index = getZorderIndex(I);
    _dataV[index] = data;
  }
  const bool GetPData(Index I) const
  {
    size_t index = getZorderIndex(I);
    return _dataP.get(index);
  }
  const void GetUData(Index I) const
  {
    size_t index = getZorderIndex(I);
    return _dataU[index];
  }
  const void GetVData(Index I) const
  {
    size_t index = getZorderIndex(I);
    return _dataV[index];
  }

  friend std::ostream& operator<<(std::ostream& o, const QuadTree<T>& tree)
  {
    o << "Printing QuadTree bitset:\n";
    o << "p:";
    o << tree.pTree.toString() << "\n";
    o << "u:";
    o << tree.uTree.toString() << "\n";
    o << "v:";
    o << tree.vTree.toString() << "\n";
    o << "data:\np:";
    o << tree._dataP.toString() << "\nu:";
    for (T* t = tree._dataU; t < tree._dataU + tree.dataAllocSize; t++)
    {
      o << (int)*t;
    }
    o << "\nv:";
    for (T* t = tree._dataV; t < tree._dataV + tree.dataAllocSize; t++)
    {
      o << (int)*t;
    }
    return o << "\n";
  }

  uint16_t hasChildrenP(size_t zorder, uint8_t d)
  {
    size_t index = (zorder >> (2 * (depth - d))) + depthOffset[d - 1];
    return pTree[index];
  }
  uint16_t hasChildrenU(size_t zorder, uint8_t d)
  {
    size_t index = (zorder >> (2 * (depth - d))) + depthOffset[d - 1];
    return uTree[index];
  }
  uint16_t hasChildrenV(size_t zorder, uint8_t d)
  {
    size_t index = (zorder >> (2 * (depth - d))) + depthOffset[d - 1];
    return vTree[index];
  }

private:
  dynamic_bitset pTree;
  dynamic_bitset uTree;
  dynamic_bitset vTree;
  size_t depth;
  dynamic_bitset _dataP;
  T* _dataU = nullptr;
  T* _dataV = nullptr;
  std::vector<size_t> depthOffset;
  size_t maxSize;
  size_t dataAllocSize;

  size_t clamp(size_t min, size_t max, size_t x)
  {
    return min > x ? min : max < x ? max
                                   : x;
  }
  bool generateGridFromImage(const std::filesystem::path& imagePath, int rgbOffset = 0)
  {
    int width, height, channels;
    unsigned char* data = stbi_load(imagePath.c_str(), &width, &height, &channels, STBI_rgb);
    assert(data);
    size_t leadingZeros = std::min(__builtin_clz(width), __builtin_clz(height)); // ggrks apparently builtinclz ignores bitwidth and is 32 bit only ty for nothing
    size_t size = ((((uint32_t)-1) >> 1) + 1) >> (leadingZeros - 1); // so 32 bit only -> cast -1 to uint32 instead of size_t
    DebugF("found image with size {}, {} and chose {}, channels {}", width, height, size, channels);
    createWithSize(size, size);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < width; y++)
      {
        int8_t down = (int8_t)data[(width * clamp(0, height - 1, y - 1) + x) * channels]; // this is boundary check down
        int8_t left = (int8_t)data[(width * y + clamp(0, width - 1, x - 1)) * channels]; // this is boundary check left
        int8_t p = (int8_t)data[(width * y + x) * channels + 0];
        int8_t u = (int8_t)data[(width * y + x) * channels + 1];
        int8_t v = (int8_t)data[(width * y + x) * channels + 2];
        if (p != 0)
        {
          SetPData({ x, y, depth });
        }
        if (u != 0)
        {
          if (left != 0)
          {
            SetUData({ x + 1, y, depth }, u);

          } else
          {
            SetUData({ x, y, depth }, u);
          }
        }

        if (v != 0)
        {
          if (down != 0)
          {
            SetVData({ x, y + 1, depth }, v);

          } else
          {
            SetVData({ x, y, depth }, v);
          }
        }
      }
    }
    // TODO: read extra u on top and extra v on right
    std::cout << "read data into tree" << std::endl;
    std::cout << *this;
    floodTreeToRoot();
    return true;
  }
  void createWithSize(size_t sizeX, size_t sizeY)
  {
    assert(isPowerOf2(sizeX) && sizeX == sizeY);
    maxSize = sizeX;
    depth = __builtin_ctz(sizeX);
    size_t allocSize = sizeX * sizeY;
    size_t tempSizeX = sizeX / 2;
    size_t tempSizeY = sizeY / 2;
    for (size_t i = 1; i < depth; ++i)
    {
      allocSize += tempSizeX * tempSizeY;
      tempSizeX /= 2;
      tempSizeY /= 2;
    }
    dataAllocSize = allocSize;
    DebugF("maxSize:{}, allocSize:{}, allocTreeSize:{}", maxSize, allocSize, allocSize - sizeX * sizeY);
    _dataP.resize(allocSize);
    _dataU = (T*)malloc(sizeof(T) * allocSize);
    _dataV = (T*)malloc(sizeof(T) * allocSize);
    std::fill(_dataU, _dataU + dataAllocSize, 0);
    std::fill(_dataV, _dataV + dataAllocSize, 0);
    pTree.resize(allocSize - sizeX * sizeY);
    uTree.resize(allocSize - sizeX * sizeY);
    vTree.resize(allocSize - sizeX * sizeY);
    depthOffset.resize(depth);
    size_t offset = 0;
    size_t nodesCount = 4;
    for (size_t i = 0; i < depth; i++)
    {
      depthOffset[i] = offset;
      offset += nodesCount;
      nodesCount *= 4;
    }
    DebugF("created tree with depth {} and offsets:{}", depth, depthOffset);
  }

#define CHILDRENAREEQUAL(dataArray, firstChild) dataArray[firstChild] == dataArray[firstChild + 1] && dataArray[firstChild + 1] == dataArray[firstChild + 2] && dataArray[firstChild + 2] == dataArray[firstChild + 3]
#define ALLCHILDRENARELEAFS(tree, firstChild) !(tree[firstChild] || tree[firstChild + 1] || tree[firstChild + 2] || tree[firstChild + 3])
  void floodTreeToRoot()
  {
    size_t offset = 1;
    size_t d = depth - 1;
    offset *= 2;
    // start by writing either data or hasChildren flag to d-1
    for (size_t y = 0; y < maxSize; y += offset)
    {
      for (size_t x = 0; x < maxSize; x += offset)
      {
        size_t indexChild = getZorderIndex({ x, y, d + 1 });
        size_t indexCurrent = getZorderIndex({ x, y, d });
        // flood P Tree
        if (CHILDRENAREEQUAL(_dataP, indexChild))
        {
          _dataP.set(indexCurrent, _dataP.get(indexChild));
        } else
        {
          pTree.set(indexCurrent, true);
        }
        // flood U Tree
        if (CHILDRENAREEQUAL(_dataV, indexChild))
        {
          _dataU[indexCurrent] = _dataU[indexChild];
        } else
        {
          uTree.set(indexCurrent, true);
        }
        // flood V Tree
        if (CHILDRENAREEQUAL(_dataV, indexChild))
        {
          _dataV[indexCurrent] = _dataV[indexChild];
        } else
        {
          vTree.set(indexCurrent, true);
        }
      }
    }

    offset = 2;
    for (d = depth - 2; d > 0; d--)
    {
      offset *= 2;
      for (size_t y = 0; y < maxSize; y += offset)
      {
        for (size_t x = 0; x < maxSize; x += offset)
        {
          size_t indexChild = getZorderIndex({ x, y, d + 1 });
          size_t indexCurrent = getZorderIndex({ x, y, d });
          // flood P tree
          if (ALLCHILDRENARELEAFS(pTree, indexChild) && CHILDRENAREEQUAL(_dataP, indexChild))
          {
            _dataP.set(indexCurrent, _dataP.get(indexChild));
          } else
          {
            pTree.set(indexCurrent, true);
          }
          // flood U Tree
          if (ALLCHILDRENARELEAFS(uTree, indexChild) && CHILDRENAREEQUAL(_dataV, indexChild))
          {
            _dataU[indexCurrent] = _dataU[indexChild];
          } else
          {
            uTree.set(indexCurrent, true);
          }
          // flood V Tree
          if (ALLCHILDRENARELEAFS(vTree, indexChild) && CHILDRENAREEQUAL(_dataV, indexChild))
          {
            _dataV[indexCurrent] = _dataV[indexChild];
          } else
          {
            vTree.set(indexCurrent, true);
          }
        }
      }
    }
  }
};

#endif // QUADTREE_H_
