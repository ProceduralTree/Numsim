#ifndef BROADCAST_H_
#define BROADCAST_H_

#include <cstdint>
#include <functional>
#include <grid/grid.h>
#include <pde/system.h>
#include <utils/index.h>
#define LOG(x) std::cout << #x << "=" << x << std::endl
inline void broadcast(
  std::function<void(PDESystem&, Index)> Operator,
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
inline void parallel_broadcast(
  std::function<void(PDESystem&, Index)> Operator,
  PDESystem& system,
  Index Begin,
  Index End)
{

#pragma omp parallel for simd collapse(2)
  for (uint16_t j = Begin.y; j <= End.y; j++)
  {
    for (uint16_t i = Begin.x; i <= End.x; i++)
    {
      Operator(system, { i, j });
    }
  }
};
constexpr void broadcast(
  std::function<void(PDESystem&, Index)> Operator,
  PDESystem& system,
  Range r)
{
  broadcast(Operator, system, r.begin, r.end);
};

constexpr void broadcast(
  std::function<void(PDESystem&, Index)> Operator,
  PDESystem& system,
  std::vector<Range> ranges)
{
  for (auto r : ranges)
    broadcast(Operator, system, r.begin, r.end);
};

constexpr void broadcast_x_boundary(
  std::function<void(PDESystem&, Index, Offset)> Operator,
  PDESystem& system, const Grid2D& grid)
{
  // TOP
  broadcast(
    [&](PDESystem& s, Index I) { Operator(s, I, Iy); },
    system, grid.end - grid.len_x, grid.end);

  // BOTTOM
  broadcast(
    [&](PDESystem& s, Index I) { Operator(s, I, -Iy); },
    system, grid.begin, grid.begin + grid.len_x);
};
constexpr void broadcast_y_boundary(
  std::function<void(PDESystem&, Index, Offset)> Operator,
  PDESystem& system, const Grid2D& grid)
{
  // LEFT
  broadcast(
    [&](PDESystem& s, Index I) { Operator(s, I, -Ix); },
    system, grid.begin, grid.begin + grid.len_y);

  // RIGHT
  broadcast(
    [&](PDESystem& s, Index I) { Operator(s, I, Ix); },
    system, grid.end - grid.len_y, grid.end);
};
constexpr void broadcast_boundary(
  std::function<void(PDESystem&, Index, Offset)> Operator,
  PDESystem& system, const Grid2D& grid)
{
  broadcast_x_boundary(Operator, system, grid);
  broadcast_y_boundary(Operator, system, grid);
};

template <typename Operator, typename... Args>
void broadcast(Operator&& O, Range r, Args&&... args)
{

#pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y; j <= r.end.y; j++)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i++)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
};

template <typename Operator, typename... Args>
void broadcast_boundary(Operator&& O, Boundaries boundaries, Args&&... args)
{
  for (auto [b, o] : boundaries)
  {
    broadcast(std::forward<Operator>(O), b, o, std::forward<Args>(args)...);
  };
};

constexpr void copy(Index I, Offset O, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };
#endif // BROADCAST_H_
