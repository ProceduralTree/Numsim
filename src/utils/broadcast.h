#ifndef BROADCAST_H_
#define BROADCAST_H_

#include <cstdint>
#include <functional>
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

// template <typename Operator, typename... Args>
//  void broadcast_blackred(Operator&& O, Range r, Args&&... args)
//{
//  #pragma omp parallel for simd collapse(2)
//    for (uint16_t j = r.begin.y; j <= r.end.y; j++)
//    {
//      for (uint16_t i = r.begin.x; i <= r.end.x; i = i + 2)
//      {
//        std::forward<Operator>(O)(Index { static_cast<uint16_t>(i + j % 2), j }, std::forward<Args>(args)...);
//      }
//    }
//  #pragma omp parallel for simd collapse(2)
//    for (uint16_t j = r.begin.y; j <= r.end.y; j++)
//    {
//      for (uint16_t i = r.begin.x; i <= r.end.x; i = i + 2)
//      {
//        std::forward<Operator>(O)(Index { static_cast<uint16_t>(i + !(j % 2)), j }, std::forward<Args>(args)...);
//      }
//    }
//  };

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
void broadcast_boundary(Operator&& O, Boundaries boundaries, Args&&... args)
{
  for (auto [b, o] : boundaries)
  {
    broadcast(std::forward<Operator>(O), b, o, std::forward<Args>(args)...);
  };
};

inline void copy(Index I, Offset O, const Grid2D& from, Grid2D& to) { to[I] = from[I]; };
#endif // BROADCAST_H_
