#include "profiler.h"
#include <chrono>
#include <fstream>
#include <stack>

#ifdef DEBUG
static std::ofstream file;

namespace Profiler {

std::stack<TimeStamp> stack;

void Init()
{
  file = std::ofstream("Profiler.csv");
}
void Push(const std::string& name)
{
  stack.emplace(name);
}
void Count()
{
  stack.top().counter++;
}
void Pop()
{
  stack.pop();
}
void Close()
{
  file.flush();
  file.close();
}
TimeStamp::TimeStamp(const std::string& name)
{
  this->name = name;
  this->start = std::chrono::high_resolution_clock::now();
}

TimeStamp::~TimeStamp()
{
  this->end = std::chrono::high_resolution_clock::now();
  file << name << ',' << start << ',' << end << ',' << end - start << ',' << counter << '\n';
}

StackHelper::StackHelper(const std::string& name)
{
  stack.emplace(name);
}
StackHelper::~StackHelper()
{
  stack.pop();
}
}
#endif
