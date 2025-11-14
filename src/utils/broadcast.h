#ifndef BROADCAST_H_
#define BROADCAST_H_

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
};

template <typename OPERATOR>
inline void broadcast(
  OPERATOR Operator,
  PDESystem& system,
  std::vector<Range> ranges)
{
  for (auto r : ranges)
    broadcast(Operator, system, r.begin, r.end);
};

template <typename Operator, typename... Args>
void broadcast_blackred(Operator&& O, Range r, Args&&... args)
{

#pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y; j <= r.end.y; j = j + 2)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i = i + 2)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
#pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y; j <= r.end.y; j = j + 2)
  {
    for (uint16_t i = r.begin.x + 1; i <= r.end.x; i = i + 2)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
#pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y + 1; j <= r.end.y; j = j + 2)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i = i + 2)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
#pragma omp parallel for simd collapse(2)
  for (uint16_t j = r.begin.y + 1; j <= r.end.y; j = j + 2)
  {
    for (uint16_t i = r.begin.x + 1; i <= r.end.x; i = i + 2)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
};

template <typename Operator, typename... Args>
void broadcast(Operator&& O, Range r, Args&&... args)
{

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

#pragma omp parallel for schedule(dynamic)
  for (uint16_t j = r.begin.y; j <= r.end.y; j++)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i++)
    {
      std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
};

template <typename Operator, typename... Args>
void broadcast_boundary(Operator&& O, Boundaries b, Args&&... args)
{
  parallel_broadcast(std::forward<Operator>(O), b.top, -Iy, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.bottom, Iy, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.left, Ix, std::forward<Args>(args)...);
  parallel_broadcast(std::forward<Operator>(O), b.right, -Ix, std::forward<Args>(args)...);
};

inline void copy(Index I, Offset O, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };
// inline void copy(Index I, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };

inline void set(Index I, Offset O, Grid2D& array, double value) { array[I] = value; };
#endif // BROADCAST_H_
