#ifndef GRID_H_
#define GRID_H_

#include <cstdint>
#include <iostream>
#include <ostream>
#include <vector>
class Grid2D {
  private:
    std::vector<double> _data;

  public:
    const uint32_t size_x;
    const uint32_t size_y;

    Grid2D(uint32_t x, uint32_t y);
    double &operator[](uint32_t x, uint32_t y);
    const double &operator[](uint32_t x, uint32_t y) const;
};
std::ostream &operator<<(std::ostream &os, const Grid2D &obj);
uint64_t interleave(uint32_t x, uint32_t y);

void laplace(const Grid2D &in, Grid2D &out);
void hello();

#endif // GRID_H_
