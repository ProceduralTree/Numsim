#ifndef BROADCAST_H_
#define BROADCAST_H_

#include "utils/partitioning.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <cstdint>
#include <grid/grid.h>
#include <pde/system.h>
#include <utils/index.h>

template <typename Operator, typename... Args>
void broadcast_blackred(Operator&& O, int parity, Range r, Args&&... args)
{
  ProfileScope("Black Iteration");
  constexpr uint16_t BLOCK_SIZE_X = 32;
  constexpr uint16_t BLOCK_SIZE_Y = 32;

  // #pragma omp loop bind(parallel) collapse(2)
  for (uint16_t by = r.begin.y; by <= r.end.y; by += BLOCK_SIZE_Y)
  {
    for (uint16_t bx = r.begin.x; bx <= r.end.x; bx += BLOCK_SIZE_X)
    {
      uint16_t y_max = std::min<uint16_t>(by + BLOCK_SIZE_Y - 1, r.end.y);
      uint16_t x_max = std::min<uint16_t>(bx + BLOCK_SIZE_X - 1, r.end.x);

      for (uint16_t j = by + parity; j <= y_max; j = j + 2)
      {
        for (uint16_t i = bx; i <= x_max; i = i + 2)
        {
          std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
        }
      }
      for (uint16_t j = by + 1 - parity; j <= y_max; j = j + 2)
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
void broadcast_red(Operator&& O, Range r, Args&&... args)
{
  ProfileScope("Red Iteration");
  constexpr uint16_t BLOCK_SIZE_X = 32;
  constexpr uint16_t BLOCK_SIZE_Y = 32;
  // #pragma omp loop bind(parallel) collapse(2)
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
  ProfileScope("Sequential Broadcast");

  // #pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y; j <= r.end.y; j++)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i++)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
};

template <typename Operator, size_t S, typename... Args>
void broadcast(Operator&& O, std::array<std::tuple<Range, Offset>, S> ranges, Args&&... args)
{
  for (auto [r, _] : ranges)
  {
    broadcast(std::forward<Operator>(O), r, std::forward<Args>(args)...);
  }
}
template <typename Operator, size_t S, typename... Args>
inline void broadcast_with_offset(Operator&& O, const std::array<std::tuple<Range, Offset>, S>& ranges, Args&&... args)
{
  for (int i = 0; i < 4; i++)
  {
    auto [r, o] = ranges[i];
    if (Settings::get().mpi.neighbours()[i][0] >= 0)
      broadcast(std::forward<Operator>(O), r - o, std::forward<Args>(args)...);
  }
}

template <typename Operator, typename... Args>
void parallel_broadcast(Operator&& O, Range r, Args&&... args)
{
  ProfileScope("Parallel Broadcast");
  constexpr uint16_t BLOCK_SIZE_X = 16;
  constexpr uint16_t BLOCK_SIZE_Y = 16;
  // #pragma omp for collapse(2)
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
  ProfileScope("Jacoby Broadcast");
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
void broadcast_boundary(Operator&& O, Partitioning::MPIInfo partitioning, Boundaries b, Args&&... args)
{
  ProfileScope("Boundary");
  if (partitioning.top_neighbor < 0)
    broadcast(std::forward<Operator>(O), b.top, -Iy, std::forward<Args>(args)...);
  if (partitioning.bottom_neighbor < 0)
    broadcast(std::forward<Operator>(O), b.bottom, Iy, std::forward<Args>(args)...);
  if (partitioning.left_neighbor < 0)
    broadcast(std::forward<Operator>(O), b.left, Ix, std::forward<Args>(args)...);
  if (partitioning.right_neighbor < 0)
    broadcast(std::forward<Operator>(O), b.right, -Ix, std::forward<Args>(args)...);
};

template <typename Operator, typename... Args>
void broadcast_ghosts(Operator&& O, Partitioning::MPIInfo partitioning, Boundaries b, Args&&... args)
{
  ProfileScope("Ghost");
  if (partitioning.top_neighbor >= 0)
    broadcast(std::forward<Operator>(O), b.top, std::forward<Args>(args)...);
  if (partitioning.bottom_neighbor >= 0)
    broadcast(std::forward<Operator>(O), b.bottom, std::forward<Args>(args)...);
  if (partitioning.left_neighbor >= 0)
    broadcast(std::forward<Operator>(O), b.left, std::forward<Args>(args)...);
  if (partitioning.right_neighbor >= 0)
    broadcast(std::forward<Operator>(O), b.right, std::forward<Args>(args)...);
};
template <typename Operator, typename... Args>
void broadcast_halo(Operator&& O, Boundaries b, Args&&... args)
{
  ProfileScope("Halo");
  broadcast(std::forward<Operator>(O), b.top, -Iy, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), b.bottom, Iy, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), b.left, Ix, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), b.right, -Ix, std::forward<Args>(args)...);
};
inline void copy(Index I, Offset O, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };
// inline void copy(Index I, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };

inline void set(Index I, Offset O, Grid2D& array, double value) { array[I] = value; };
#endif // BROADCAST_H_
