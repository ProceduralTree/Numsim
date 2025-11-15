#ifndef PROFILER_H_
#define PROFILER_H_

#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>

struct Timing
{
  std::string name;
  int calls = 0;
  double totalTime = 0;

  void add(double dt);
  double averageTime() const;
};

class Profiler
{
public:
  void addTiming(const std::string& name, long long dt);
  void print() const;

  static Profiler& instance(); // singleton pattern for global access

private:
  mutable std::mutex mutex_;
  std::unordered_map<std::string, Timing> timingsMap_;
};

class Scope
{
public:
  explicit Scope(const std::string& name);
  ~Scope();

private:
  std::string name_;
  std::chrono::high_resolution_clock::time_point start_;
};

#endif // PROFILER_H_
