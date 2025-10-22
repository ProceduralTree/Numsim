#include "indexing.h"
#include <cstdint>
#include <immintrin.h>
#include <tuple>

uint32_t top(uint32_t z) {
    uint32_t tmp = z & maskx;
    tmp--;
    tmp &= maskx;
    tmp |= z & masky;
    return tmp;
}
uint32_t bottom(uint32_t z) {
    uint32_t tmp = z & masky;
    tmp++;
    tmp &= maskx;
    tmp |= z & masky;
    return tmp;
}
uint32_t left(uint32_t z) {
    uint32_t tmp = z & masky;
    tmp--;
    tmp &= masky;
    tmp |= z & maskx;
    return tmp;
}
uint32_t right(uint32_t z) {
    uint32_t tmp = z & maskx;
    tmp++;
    tmp &= masky;
    tmp |= z & maskx;
    return tmp;
}

uint32_t interleave(uint16_t x, uint16_t y) {
    uint32_t masky = 0xAAAA;
    uint32_t maskx = 0x5555;
    uint32_t xbits = _pdep_u32(x, maskx);
    uint32_t ybits = _pdep_u32(y, masky);
    return xbits | ybits;
}

std::tuple<uint32_t, uint32_t> detangle(uint32_t zindex) {
    uint32_t masky = 0xAAAA;
    uint32_t maskx = 0x5555;
    uint32_t x = _pext_u32(zindex, maskx);
    uint32_t y = _pext_u32(zindex, masky);
    return {x, y};
}
