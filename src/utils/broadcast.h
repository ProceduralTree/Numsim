#ifndef BROADCAST_H_
#define BROADCAST_H_

#include "utils/profiler.h"
#include <cstdint>
#include <grid/grid.h>
#include <pde/system.h>
#include <utils/index.h>
#define LOG(x) std::cout << #x << "=" << x << std::endl
template <typename OPERATOR>
inline void broadcast(
  OPERATOR Operator,
  PDESystem& system,
  Index Begin,
  Index End)
{

  // #pragma omp parallel for simd collapse(2)
  for (uint16_t j = Begin.y; j <= End.y; j++)
  {
    for (uint16_t i = Begin.x; i <= End.x; i++)
    {
      Operator(system, { i, j });
    }
  }
};
template <typename OPERATOR>
inline void broadcast(
  OPERATOR Operator,
  PDESystem& system,
  Range r)
{
  broadcast(Operator, system, r.begin, r.end);
}

template <typename OPERATOR>
inline void broadcast(
  OPERATOR Operator,
  PDESystem& system,
  std::vector<Range> ranges)
{
  for (auto r : ranges)
    broadcast(Operator, system, r.begin, r.end);
}

template <typename Operator, typename... Args>
void broadcast_blackred(Operator&& O, Range r, Args&&... args)
{
  Scope scope("Black Red Iteration");
  constexpr uint16_t BLOCK_SIZE_X = 32;
  constexpr uint16_t BLOCK_SIZE_Y = 32;

#pragma omp loop bind(parallel) collapse(2)
  for (uint16_t by = r.begin.y; by <= r.end.y; by += BLOCK_SIZE_Y)
  {
    for (uint16_t bx = r.begin.x; bx <= r.end.x; bx += BLOCK_SIZE_X)
    {
      uint16_t y_max = std::min<uint16_t>(by + BLOCK_SIZE_Y - 1, r.end.y);
      uint16_t x_max = std::min<uint16_t>(bx + BLOCK_SIZE_X - 1, r.end.x);

      for (uint16_t j = by; j <= y_max; j = j + 2)
      {
        for (uint16_t i = bx; i <= x_max; i = i + 2)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
      for (uint16_t j = by + 1; j <= y_max; j = j + 2)
      {
        for (uint16_t i = bx + 1; i <= x_max; i = i + 2)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
    }
  }
#pragma omp loop bind(parallel) collapse(2)
  for (uint16_t by = r.begin.y; by <= r.end.y; by += BLOCK_SIZE_Y)
  {
    for (uint16_t bx = r.begin.x; bx <= r.end.x; bx += BLOCK_SIZE_X)
    {
      uint16_t y_max = std::min<uint16_t>(by + BLOCK_SIZE_Y - 1, r.end.y);
      uint16_t x_max = std::min<uint16_t>(bx + BLOCK_SIZE_X - 1, r.end.x);

      for (uint16_t j = by + 1; j <= y_max; j = j + 2)
      {
        for (uint16_t i = bx; i <= x_max; i = i + 2)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
      for (uint16_t j = by; j <= y_max; j = j + 2)
      {
        for (uint16_t i = bx + 1; i <= x_max; i = i + 2)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
    }
  }
}

template <typename Operator, typename... Args>
void broadcast(Operator&& O, Range r, Args&&... args)
{
  ProfileScope("Broadcast");

  Scope scope("Sequential Broadcast");
  // #pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y; j <= r.end.y; j++)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i++)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
};
template <typename Operator, typename... Args>
void parallel_broadcast(Operator&& O, Range r, Args&&... args)
{
  Scope scope("Parallel Broadcast");
  constexpr uint16_t BLOCK_SIZE_X = 16;
  constexpr uint16_t BLOCK_SIZE_Y = 16;
#pragma omp parallel loop bind(parallel) collapse(2)
  for (uint16_t by = r.begin.y; by <= r.end.y; by += BLOCK_SIZE_Y)
  {
    for (uint16_t bx = r.begin.x; bx <= r.end.x; bx += BLOCK_SIZE_X)
    {
      uint16_t y_max = std::min<uint16_t>(by + BLOCK_SIZE_Y - 1, r.end.y);
      uint16_t x_max = std::min<uint16_t>(bx + BLOCK_SIZE_X - 1, r.end.x);
      for (uint16_t j = by; j <= y_max; j++)
      {
#pragma omp simd
        for (uint16_t i = bx; i <= x_max; i++)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
    }
  }
};

template <typename Operator, typename... Args>
void test_broadcast(Operator&& O, Range r, Args&&... args)
{
  Scope scope("Jacoby Broadcast");
  constexpr uint16_t BLOCK_SIZE_X = 16;
  constexpr uint16_t BLOCK_SIZE_Y = 16;
  // #pragma omp parallel for schedule(static) collapse(2)
  for (uint16_t by = r.begin.y; by <= r.end.y; by += BLOCK_SIZE_Y)
  {
    for (uint16_t bx = r.begin.x; bx <= r.end.x; bx += BLOCK_SIZE_X)
    {
      uint16_t y_max = std::min<uint16_t>(by + BLOCK_SIZE_Y - 1, r.end.y);
      uint16_t x_max = std::min<uint16_t>(bx + BLOCK_SIZE_X - 1, r.end.x);
      for (uint16_t j = by; j <= y_max; j++)
      {
#pragma omp simd
        for (uint16_t i = bx; i <= x_max; i++)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
    }
  }
};
template <typename Operator, typename... Args>
void broadcast_boundary(Operator&& O, Boundaries b, Args&&... args)
{
  Scope scope("Boundary");
  parallel_broadcast(std::forward<Operator>(O), b.top, -Iy, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.bottom, Iy, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.left, Ix, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.right, -Ix, std::forward<Args>(args)...);
};

inline void copy(Index I, Offset O, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };
// inline void copy(Index I, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };

inline void set(Index I, Offset O, Grid2D& array, double value) { array[I] = value; };
#endif // BROADCAST_H_
