#ifndef GRID_H_
#define GRID_H_

#include <cstdint>
#include <iostream>
#include <ostream>
#include <vector>

#include "indexing.h"

template <Indexing I> class Grid2D {
  private:
    std::vector<double> _data;

  public:
    const uint16_t size_x;
    const uint16_t size_y;

    Grid2D(uint16_t x, uint16_t y);

    double &operator[](uint16_t x, uint16_t y);
    const double &operator[](uint16_t x, uint16_t y) const;
    double &operator[](uint32_t z);
    const double &operator[](uint32_t z) const;
    double &operator()(uint16_t x, uint16_t y);
    const double &operator()(uint16_t x, uint16_t y) const;

    const uint32_t elements() const { return this->_data.size(); }
};
template <Indexing I> void laplace(const Grid2D<I> &in, Grid2D<I> &out);

template <>
void laplace(const Grid2D<Indexing::Cartesian> &in,
             Grid2D<Indexing::Cartesian> &out);
template <>
void laplace(const Grid2D<Indexing::ZOrder> &in, Grid2D<Indexing::ZOrder> &out);

template <Indexing I>
std::ostream &operator<<(std::ostream &os, const Grid2D<I> &obj);

template class Grid2D<Indexing::Cartesian>;
template class Grid2D<Indexing::ZOrder>;

#endif // GRID_H_
