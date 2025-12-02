#ifndef DISTRIBUTED_H_
#define DISTRIBUTED_H_
#include "grid/grid.h"
#include "utils/Logger.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <grid/grid.h>
#include <mpi.h>
#include <strings.h>
#include <utility>
#include <utils/broadcast.h>
#include <utils/partitioning.h>

inline size_t len(Range r)
{
  return (r.end.x - r.begin.x + 1) * (r.end.y - r.begin.y + 1);
};

struct MPI_COMM_BUFFER
{
  MPI_Comm comm;
  Grid2D& comm_array;
  std::array<std::tuple<Range, Offset>, 4> communication_boundary;
  std::array<MPI_Request, 4> request;
  std::array<double*, 4> sendbuffer;
  std::array<double*, 4> recivebuffer;

  MPI_COMM_BUFFER(Grid2D& comm_array, std::array<std::tuple<Range, Offset>, 4> ghosts, MPI_Comm comm, MPIInfo& info)
    : comm(comm)
    , comm_array(comm_array)
    , communication_boundary(ghosts)
    , request({ MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL })
    , sendbuffer()
    , recivebuffer()
  {

    for (int i = 0; i < 4; i++)
    {
      // DebugF("Neighbour: {} , Rank {}", info.neighbours()[i], info.rank);

      if (info.neighbours()[i][0] >= 0)
      {

        ProfileScope("MPI Communication Init");

        auto [r, o] = ghosts[i];
        // DebugF("Allocating sendrevieve buffers of size {} from rank {}", len(r), info.rank);
        // DebugF("Range : {{x={} , y={} }} -> {{x={} , y={} }}", r.begin.x, r.begin.y, r.end.x, r.end.y);
        // DebugF("Offset : {{x={} , y={} }}", o.x, o.y);
        sendbuffer[i] = (double*)malloc(len(r) * sizeof(double));
        recivebuffer[i] = (double*)malloc(len(r) * sizeof(double));
        comm_array.get(sendbuffer[i], r - o);
        // DebugF("Request {}", request[i]);
        MPI_Isendrecv(
          sendbuffer[i],
          len(r),
          MPI_DOUBLE,
          info.neighbours()[i][0],
          i,
          recivebuffer[i],
          len(r),
          MPI_DOUBLE,
          info.neighbours()[i][0],
          info.neighbours()[i][1],
          comm,
          &request[i]);
      }
    }
  }
  ~MPI_COMM_BUFFER()
  {
    ProfileScope("MPI Communication Wait");

    // guaranties a maximum of 4 iterations(while is evil)
    for (int i = 0; i < 4; i++)
    {
      int indices[4];
      int outcout;

      // scipped if all requests are MPI_REQUEST_NULL with outcout=MPI_UNDEFINED
      MPI_Waitsome(4, request.data(), &outcout, indices, MPI_STATUSES_IGNORE);
      if (outcout == MPI_UNDEFINED)
        break;
      // scipped if outcount=0
      for (int succes = 0; succes < outcout; succes++)
      {
        int index = indices[succes];
        // request[index] = MPI_REQUEST_NULL;
        auto [r, o] = communication_boundary[index];
        comm_array.set(recivebuffer[index], r);
        free(recivebuffer[index]);
        free(sendbuffer[index]);
      }
    }
  }
};

template <typename Operator, typename... Args>
void distributed_broadcast(Operator&& O, MPIInfo p, Range r, Grid2D& comm_array, Args&&... args)
{
  assert(r.end.x - r.begin.x > 2);
  assert(r.end.y - r.begin.y > 2);
  Range inner = Range { r.begin + II, r.end - II };
  Boundaries border = Boundaries(inner.begin, inner.end);
  Boundaries ghosts = Boundaries(r.begin, r.end);

  //  copy boundary sendbuff
  // broadcast(std::forward<Operator>(O), r, std::forward<Args>(args)...);
  // broadcast(std::forward<Operator>(O), inner, std::forward<Args>(args)...);
  broadcast(std::forward<Operator>(O), ghosts, std::forward<Args>(args)...);
  MPI_COMM_BUFFER* comm_buffer = new MPI_COMM_BUFFER(comm_array, ghosts.all, MPI_COMM_WORLD, p);
  broadcast(std::forward<Operator>(O), r, std::forward<Args>(args)...);
  delete comm_buffer;
};

#endif // DISTRIBUTED_H_
