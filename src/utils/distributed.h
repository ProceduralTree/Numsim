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
  const Partitioning::MPIInfo& info;
  MPI_Comm comm;
  Grid2D& comm_array;
  std::array<std::tuple<Range, Offset>, 4> communication_boundary;
  std::array<MPI_Request, 4> requestS;
  std::array<MPI_Request, 4> requestR;
  std::array<double*, 4> sendbuffer;
  std::array<double*, 4> recivebuffer;

  MPI_COMM_BUFFER(Grid2D& comm_array, std::array<std::tuple<Range, Offset>, 4> ghosts, MPI_Comm comm, Partitioning::MPIInfo& info, int id = 0)
    : info(info)
    , comm(comm)
    , comm_array(comm_array)
    , communication_boundary(ghosts)
    , requestS({ MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL })
    , requestR({ MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL })
    , sendbuffer()
    , recivebuffer()
  {
    Init(ghosts);
  };
  ~MPI_COMM_BUFFER()
  {
    // not needed
    for (size_t i = 0; i < 4; i++)
    {
      free(sendbuffer[i]);
      free(recivebuffer[i]);
    }
  }
  void Init(const std::array<std::tuple<Range, Offset>, 4>& ghosts)
  {
    for (int i = 0; i < 4; i++)
    {
      if (info.neighbours()[i][0] >= 0)
      {
        auto [r, o] = ghosts[i];
        recivebuffer[i] = (double*)malloc(len(r) * sizeof(double));
      }
    }

    for (int i = 0; i < 4; i++)
    {
      if (info.neighbours()[i][0] >= 0)
      {
        auto [r, o] = ghosts[i];
        sendbuffer[i] = (double*)malloc(len(r) * sizeof(double));
      }
    }
  }
  void Send(const std::array<std::tuple<Range, Offset>, 4>& ghosts, int id = 0)
  {
    ProfileScope("MPI Communication Init");
    for (int i = 0; i < 4; i++)
    {
      if (info.neighbours()[i][0] >= 0)
      {
        auto [r, o] = ghosts[i];
        MPI_Irecv(recivebuffer[i], len(r), MPI_DOUBLE, info.neighbours()[i][0], info.neighbours()[i][1] + id, comm, &requestR[i]);
      }
    }

    for (int i = 0; i < 4; i++)
    {
      if (info.neighbours()[i][0] >= 0)
      {
        auto [r, o] = ghosts[i];
        comm_array.get(sendbuffer[i], r - o);
        MPI_Isend(sendbuffer[i], len(r), MPI_DOUBLE, info.neighbours()[i][0], i + id, comm, &requestS[i]);
      }
    }
  }
  void Receive()
  {
    ProfileScope("MPI Communication Wait");

    // guaranties a maximum of 4 iterations(while is evil)
    for (int i = 0; i < 4; i++)
    {
      int indices[4];
      int outcout;

      // scipped if all requests are MPI_REQUEST_NULL with outcout=MPI_UNDEFINED
      MPI_Waitsome(4, requestR.data(), &outcout, indices, MPI_STATUSES_IGNORE);
      if (outcout == MPI_UNDEFINED)
        break;
      // scipped if outcount=0
      for (int succes = 0; succes < outcout; succes++)
      {
        int index = indices[succes];
        // request[index] = MPI_REQUEST_NULL;
        auto [r, o] = communication_boundary[index];
        comm_array.set(recivebuffer[index], r);
        // free(recivebuffer[index]);
      }
    }
    for (int i = 0; i < 4; i++)
    {
      int indices[4];
      int outcout;

      // scipped if all requests are MPI_REQUEST_NULL with outcout=MPI_UNDEFINED
      MPI_Waitsome(4, requestS.data(), &outcout, indices, MPI_STATUSES_IGNORE);
      if (outcout == MPI_UNDEFINED)
        break;
      // scipped if outcount=0
      for (int succes = 0; succes < outcout; succes++)
      {
        int index = indices[succes];
        // free(sendbuffer[index]);
      }
    }
  }
};

template <typename Operator, typename... Args>
void distributed_broadcast(Operator&& O, Partitioning::MPIInfo p, Range r, Grid2D& comm_array, Args&&... args)
{
  assert(r.end.x - r.begin.x > 2);
  assert(r.end.y - r.begin.y > 2);
  Range inner = Range { r.begin + II, r.end - II };
  Boundaries border = Boundaries(inner.begin, inner.end);
  Boundaries ghosts = Boundaries(r.begin, r.end);

  //  copy boundary sendbuff
  broadcast(std::forward<Operator>(O), border.unique(), std::forward<Args>(args)...);
  // broadcast(std::forward<Operator>(O), r, std::forward<Args>(args)...);
  MPI_COMM_BUFFER comm_buffer(comm_array, ghosts.all, MPI_COMM_WORLD, p);
  comm_buffer.Send(ghosts.all);
  broadcast(std::forward<Operator>(O), inner, std::forward<Args>(args)...);
  comm_buffer.Receive();
};

#endif // DISTRIBUTED_H_
