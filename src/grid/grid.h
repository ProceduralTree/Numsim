#ifndef GRID_H_
#define GRID_H_

#include <cstdint>
#include <iostream>
#include <ostream>
#include <vector>

#include "indexing.h"

#define CARTESIAN

class Grid2D
{
private:
  std::vector<double> _data;

public:
  const uint16_t size_x;
  const uint16_t size_y;

  Grid2D(uint16_t x, uint16_t y);

  double& operator[](uint16_t x, uint16_t y);
  const double& operator[](uint16_t x, uint16_t y) const;
  double& operator[](uint32_t z);
  const double& operator[](uint32_t z) const;

  const uint32_t elements() const { return this->_data.size(); }
};

void laplace(const Grid2D& in, Grid2D& out);
std::ostream& operator<<(std::ostream& os, const Grid2D& obj);

#endif // GRID_H_
