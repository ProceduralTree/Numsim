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

  size_t getZorderIndexWithoutDepthCorrection(const Index I, const std::vector<size_t>& depthOffset) const
  {
    size_t zorder = IndexToZOrder_l(I.x, I.y);
    return (zorder) + depthOffset[I.depth - 1];
  }
  std::string toString(const size_t depth, const std::vector<size_t>& depthOffset) const
  {
    std::stringstream ss;
    size_t size = 2;
    for (size_t d = 1; d <= depth; d++)
    {
      ss << d << ":\n";
      for (size_t y = 0; y < size; y++)
      {
        for (size_t x = 0; x < size; x++)
        {
          ss << operator[](getZorderIndexWithoutDepthCorrection({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d) }, depthOffset));
        }
        ss << "\n";
      }
      size *= 2;
    }
    return ss.str();
  }

private:
  uint8_t* _data = nullptr;
  size_t _numBits = 0;
  size_t _numBytes = 0;
};

struct QuadTree
{
  // for image import the following values will be converted to the following boundary types:
  // 0 = OUTSIDE 1,2,3,4,5 = INSIDE,BOT,TOP,LEFT,RIGHT
  //
  enum class CellType : uint8_t
  {
    // sides
    OUTSIDE = 0b0000'0000,
    INSIDE = 0b0000'0001,
    BOTTOM = 0b0000'0010,
    TOP = 0b0000'0100,
    LEFT = 0b0000'1000,
    RIGHT = 0b0001'0000,
    MIXED = 0b0010'0000, // used for building tree, if we merged 4 cells with different types we will set this flag
    BOUNDARYMASK = INSIDE | BOTTOM | TOP | LEFT | RIGHT
  };
  QuadTree(const std::filesystem::path& imagePath)
  {
    generateGridFromImage(imagePath);
  }
  QuadTree(size_t sizeX, size_t sizeY)
  {
    createWithSize(sizeX, sizeY);
  }
  ~QuadTree()
  {
    free(_dataU);
    _dataU = nullptr;
    free(_dataV);
    _dataV = nullptr;
  }
  QuadTree(const QuadTree&) = delete;
  QuadTree& operator=(const QuadTree&) = delete;
  QuadTree(QuadTree&&) = delete;
  QuadTree& operator=(QuadTree&&) = delete;

  size_t getZorderIndex(Index I) const
  {
    size_t zorder = IndexToZOrder_l(I.x, I.y);
    return (zorder >> (2 * (depth - I.depth))) + depthOffset[I.depth - 1];
  }
  size_t getZorderIndexWithoutDepthCorrection(const Index I, const std::vector<size_t>& depthOffset) const
  {
    size_t zorder = IndexToZOrder_l(I.x, I.y);
    return (zorder) + depthOffset[I.depth - 1];
  }
  const uint8_t& getULowestValue(Index I) const
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
  void SetPData(Index I, CellType flag)
  {
    size_t index = getZorderIndex(I);
    _dataP[index] = flag;
  }
  void SetUData(Index I, uint8_t data)
  {
    size_t index = getZorderIndex(I);
    _dataU[index] = data;
  }
  void SetVData(Index I, uint8_t data)
  {
    size_t index = getZorderIndex(I);
    _dataV[index] = data;
  }
  const CellType GetPData(Index I) const
  {
    if (I.depth == 0)
    {
      return (CellType)((uint8_t)_dataP[0] | (uint8_t)_dataP[1] | (uint8_t)_dataP[2] | (uint8_t)_dataP[3]);
    }
    size_t index = getZorderIndex(I);
    return _dataP[index];
  }
  const uint8_t GetUData(Index I) const
  {
    size_t index = getZorderIndex(I);
    return _dataU[index];
  }
  const uint8_t GetVData(Index I) const
  {
    size_t index = getZorderIndex(I);
    return _dataV[index];
  }

  friend std::ostream& operator<<(std::ostream& o, const QuadTree& tree)
  {
    o << "Printing QuadTree bitset:\n";
    o << "p:";
    o << tree.pTree.toString(tree.depth - 1, tree.depthOffset) << "\n";
    // o << "u:";
    // o << tree.uTree.toString(tree.depth - 1, tree.depthOffset) << "\n";
    // o << "v:";
    // o << tree.vTree.toString(tree.depth - 1, tree.depthOffset) << "\n";
    o << "\np:";
    size_t offsetIndex = 0;
    size_t offsetCounter = 0;
    for (uint8_t* t = (uint8_t*)tree._dataP; t < (uint8_t*)tree._dataP + tree.dataAllocSize; t++, offsetCounter++)
    {
      if (tree.depthOffset[offsetIndex] == offsetCounter)
      {
        o << "\n";
        offsetIndex++;
      }
      o << (int)*t << "|";
      if (offsetCounter % 4 == 0)
        o << " ";
    }
    o << "\nu:";
    // for (uint8_t* t = tree._dataU; t < tree._dataU + tree.dataAllocSize; t++)
    // {
    //   o << (int)*t;
    // }
    // o << "\nv:";
    // for (uint8_t* t = tree._dataV; t < tree._dataV + tree.dataAllocSize; t++)
    // {
    //   o << (int)*t;
    // }
    return o << "\n";
  }

  uint16_t hasChildrenP(size_t zorder, uint8_t d)
  {
    if (d == 0)
      return true;
    if (d > depth)
      return false;
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
  size_t getDepth() const { return depth; }

private:
  dynamic_bitset pTree;
  dynamic_bitset uTree;
  dynamic_bitset vTree;
  size_t depth;
  CellType* _dataP = nullptr;
  uint8_t* _dataU = nullptr;
  uint8_t* _dataV = nullptr;
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
    size_t size = ((((uint32_t)-1) >> 1) + 1) >> (leadingZeros); // so 32 bit only -> cast -1 to uint32 instead of size_t
    DebugF("found image with size {}, {} and chose {}, channels {}", width, height, size, channels);
    createWithSize(size, size);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < width; y++)
      {
        // int8_t down = (int8_t)data[(width * clamp(0, height - 1, y - 1) + x) * channels]; // this is boundary check down
        // int8_t left = (int8_t)data[(width * y + clamp(0, width - 1, x - 1)) * channels]; // this is boundary check left
        uint8_t p = (uint8_t)data[(width * y + x) * channels + 0];
        int8_t u = (int8_t)data[(width * y + x) * channels + 1];
        int8_t v = (int8_t)data[(width * y + x) * channels + 2];
        Index I(x, y, depth);
        // 0 = OUTSIDE 1,2,3,4,5 = INSIDE,BOT,TOP,LEFT,RIGHT
        if (p != 0)
        {
          SetPData(I, (CellType)(1 << (p - 1)));
        }
        if (u != 0)
        {
          SetUData(I, u);
        }
        if (v != 0)
        {
          SetVData(I, v);
        }
      }
    }
    floodTreeToRoot();
    droughtLeafsByOne();
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
    _dataP = (CellType*)malloc(sizeof(CellType) * allocSize);
    _dataU = (uint8_t*)malloc(sizeof(uint8_t) * allocSize);
    _dataV = (uint8_t*)malloc(sizeof(uint8_t) * allocSize);
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
#define NOCHILDISMIXED(tree, firstChild) !((((uint8_t)tree[firstChild] | (uint8_t)tree[firstChild + 1] | (uint8_t)tree[firstChild + 2] | (uint8_t)tree[firstChild + 3]) & (uint8_t)CellType::MIXED) > 0)
#define MIXCHILDREN(tree, firstChild) ((uint8_t)tree[firstChild] | (uint8_t)tree[firstChild + 1] | (uint8_t)tree[firstChild + 2] | (uint8_t)tree[firstChild + 3])
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
        size_t indexChild = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d + 1) });
        size_t indexCurrent = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d) });
        // flood P Tree
        if (CHILDRENAREEQUAL(_dataP, indexChild))
        {
          _dataP[indexCurrent] = _dataP[indexChild];
        } else
        {
          pTree.set(indexCurrent, true);
          uint8_t currentP = MIXCHILDREN(_dataP, indexChild);
          _dataP[indexCurrent] = (CellType)(currentP | (uint8_t)CellType::MIXED);
        }
        // flood U Tree
        if (CHILDRENAREEQUAL(_dataU, indexChild))
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
          size_t indexChild = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d + 1) });
          size_t indexCurrent = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d) });
          // flood P tree
          if (ALLCHILDRENARELEAFS(pTree, indexChild) && CHILDRENAREEQUAL(_dataP, indexChild) && NOCHILDISMIXED(_dataP, indexChild))
          {
            _dataP[indexCurrent] = _dataP[indexChild];
          } else
          {
            pTree.set(indexCurrent, true);
            uint8_t currentP = MIXCHILDREN(_dataP, indexChild);
            _dataP[indexCurrent] = (CellType)(currentP | (uint8_t)CellType::MIXED);
          }
          // flood U Tree
          if (ALLCHILDRENARELEAFS(uTree, indexChild) && CHILDRENAREEQUAL(_dataU, indexChild))
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
  void droughtLeafsByOne()
  {
    size_t offset = 1;
    for (size_t d = depth - 1; d > 1; d--)
    {
      for (size_t y = 0; y < maxSize; y += offset)
      {
        for (size_t x = 0; x < maxSize; x += offset)
        {
          size_t index = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d) });
          size_t parentIndex = getZorderIndex({ static_cast<uint16_t>(x), static_cast<uint16_t>(y), static_cast<uint16_t>(d - 1) });
          if (pTree.get(parentIndex))
          {
            pTree.set(index);
          }
        }
      }
      offset *= 2;
    }
    pTree.set(0);
    pTree.set(1);
    pTree.set(2);
    pTree.set(3);
  }
};

#endif // QUADTREE_H_
