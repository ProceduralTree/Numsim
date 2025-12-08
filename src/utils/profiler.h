#pragma once

#define DEBUG
#include <string>

#include <chrono>
namespace Profiler {
enum Type
{
  FILE,
  ACCUMULATE,
  ACCUMULATEPAR
};
#ifdef DEBUG
struct TimeStamp
{
  std::string name;
  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point end;
  uint32_t counter;

  TimeStamp(const std::string& name);
  ~TimeStamp();

  TimeStamp(const TimeStamp&) = delete;
  TimeStamp& operator=(const TimeStamp&) = delete;
  TimeStamp(TimeStamp&& other) noexcept
  {
    name = std::move(other.name);
    start = std::move(other.start);
    end = std::move(other.end);
    counter = other.counter;
  }
  TimeStamp& operator=(TimeStamp&& other) noexcept
  {
    name = std::move(other.name);
    start = std::move(other.start);
    end = std::move(other.end);
    counter = other.counter;
    return *this;
  }

  void Begin();
  void End();
};
struct StackHelper
{
  uint32_t number;
  StackHelper(const std::string& name);
  ~StackHelper();
};
#endif

#ifdef DEBUG
void Init(Type type);
void Push(const std::string& name);
void Count();
void Pop();
void Close();
void PrintStack();
#else
inline void Init(Type type) { }
inline void Push(const std::string& name) { }
inline void Pop() { }
inline void Close() { }
inline void PrintStack() { }
#endif

#ifdef DEBUG
#define CONCAT(a, b) a##b
#define ProfileScopeI(name, line) Profiler::StackHelper CONCAT(stackHelper, line) = Profiler::StackHelper(name)
#define ProfileScope(name) ProfileScopeI(name, __LINE__)
#define ProfileCount() Profiler::Count()
#define ProfilePush(name) Profiler::Push(name);
#define ProfilePop() Profiler::Pop();
#else
#define ProfileScope(name)
#define ProfileCount()
#define ProfilePush(name)
#define ProfilePop()
#endif

}
#undef DEBUG
