#include "profiler.h"
#include "Logger.h"
#include "settings.h"
#include <chrono>
#include <fstream>
#include <stack>
#include <unordered_map>

#ifdef DEBUG
static std::ofstream file;
static Profiler::Type profType;

struct Accumulated
{
  size_t count;
  std::chrono::high_resolution_clock::duration min;
  std::chrono::high_resolution_clock::duration duration;
  std::chrono::high_resolution_clock::duration max;
};
static std::unordered_map<std::string, Accumulated> map;

namespace Profiler {

std::stack<TimeStamp> stack;

void Init(Type type)
{
  profType = type;
  switch (profType)
  {
  case FILE:
  default:
    file = std::ofstream("Profiler.csv");
    break;
  case ACCUMULATE:
    break;
  case ACCUMULATEPAR:
    break;
  }
}
void Push(const std::string& name)
{
  stack.emplace(name);
  stack.top().Begin();
}
void Count()
{
  stack.top().counter++;
}
void Pop()
{
  stack.top().End();
  stack.pop();
}
void Close()
{
  switch (profType)
  {
  case FILE:
  default:
    file.flush();
    file.close();
    break;
  case ACCUMULATE:
    file = std::ofstream("Profiler.csv");

    for (const auto& m : map)
    {
      const auto& s = m.first;
      const auto& v = m.second;

      file << s << ',' << v.count << ',' << v.duration << ',' << v.duration / v.count << ',' << v.min << ',' << v.max << '\n';
    }

    file.flush();
    file.close();
    break;
  case ACCUMULATEPAR:
    file = std::ofstream(std::format("Profiler{}.csv", Settings::get().mpi.rank));

    for (const auto& m : map)
    {
      const auto& s = m.first;
      const auto& v = m.second;

      file << s << ',' << v.count << ',' << v.duration << ',' << v.duration / v.count << ',' << v.min << ',' << v.max << '\n';
    }

    file.flush();
    file.close();
    //
    // create comp data from map

    // mpi gather for comp data
    // write gather data on main rank
    break;
  }
}
void PrintStack()
{
  std::stack<TimeStamp> tempStack;
  while (!stack.empty())
  {
    auto temp = std::move(stack.top());
    stack.pop();
    DebugF("{}, {}, {}", temp.name, temp.start, temp.counter);
    tempStack.push(std::move(temp));
  }
  while (!tempStack.empty())
  {
    auto temp = std::move(tempStack.top());
    tempStack.pop();
    stack.push(std::move(temp));
  }
}
TimeStamp::TimeStamp(const std::string& name)
{
  this->name = name;
}

TimeStamp::~TimeStamp()
{
}

void TimeStamp::Begin()
{
  this->start = std::chrono::high_resolution_clock::now();
}

void TimeStamp::End()
{
  this->end = std::chrono::high_resolution_clock::now();
  switch (profType)
  {
  case FILE:
  default:
    file << name << ',' << start << ',' << end << ',' << end - start << ',' << counter << '\n';
    break;
  case ACCUMULATE:
  {
    auto v = map.find(name);
    if (v == map.end())
    {
      auto currentDuration = end - start;
      map.insert({ name, { 1, currentDuration, currentDuration, currentDuration } });
    } else
    {
      v->second.count++;
      auto currentDuration = end - start;
      v->second.duration += currentDuration;
      v->second.min = std::min(v->second.min, currentDuration);
      v->second.max = std::max(v->second.max, currentDuration);
    }
  }
  break;
  case ACCUMULATEPAR:
  {
    auto v = map.find(name);
    if (v == map.end())
    {
      auto currentDuration = end - start;
      map.insert({ name, { 1, currentDuration, currentDuration, currentDuration } });
    } else
    {
      v->second.count++;
      auto currentDuration = end - start;
      v->second.duration += currentDuration;
      v->second.min = std::min(v->second.min, currentDuration);
      v->second.max = std::max(v->second.max, currentDuration);
    }
  }
  break;
  }
}

StackHelper::StackHelper(const std::string& name)
{
  Push(name);
}
StackHelper::~StackHelper()
{
  Pop();
}
}
#endif
