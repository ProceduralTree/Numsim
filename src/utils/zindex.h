#ifndef ZINDEX_H_
#define ZINDEX_H_
#include <immintrin.h>
#include <utils/index.h>
#include <valarray>

// named after most significant bit
#define MASK1 0xAAAA
#define MASK0 0x5555

struct ZIndex
{
  uint32_t z;

  ZIndex operator+(ZIndex I) const
  {

    uint32_t idx = (((z | MASK1) + (I.z & MASK0)) & MASK0) | (((z | MASK0) + (I.z & MASK1)) & MASK1);
    return ZIndex { idx };
  };
  ZIndex operator-(ZIndex I) const
  {
    uint32_t idx = (((z | MASK1) - (I.z & MASK0)) & MASK0) | (((z | MASK0) - (I.z & MASK1)) & MASK1);
    return ZIndex { idx };
  };
  ZIndex(Index I)
  {
    uint16_t x = _pdep_u32(I.x, MASK1);
    uint16_t y = _pdep_u32(I.y, MASK0);
    z = x | y;
  };
};

Index::Index(ZIndex I)
{
  x = _pext_u32(I.z, MASK1);
  y = _pext_u32(I.z, MASK0);
};

#endif // ZINDEX_H_
