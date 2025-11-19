#include "profiler.h"
#include "Logger.h"
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
  file.flush();
  file.close();
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
