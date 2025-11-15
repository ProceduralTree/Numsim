#include <iomanip>
#include <iostream>
#include <ratio>
#include <utils/profiler.h>

// ---------------- Timing ----------------
void Timing::add(double dt)
{
  totalTime += dt;
  calls++;
}

double Timing::averageTime() const
{
  return calls > 0 ? totalTime / calls : 0.0;
}

// ---------------- Profiler ----------------
Profiler& Profiler::instance()
{
  static Profiler profilerInstance;
  return profilerInstance;
}

void Profiler::addTiming(const std::string& name, long long dt)
{
  std::lock_guard<std::mutex> lock(mutex_);
  timingsMap_[name].name = name;
  timingsMap_[name].add(dt);
}

void Profiler::print() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  std::cout << "\n===== PROFILER RESULTS =====\n"
            << std::setprecision(5);
  for (const auto& [name, timing] : timingsMap_)
  {
    std::cout << std::setw(20) << name
              << " | calls:\t " << std::setw(10) << timing.calls
              << " | total:\t " << timing.totalTime / 1e6 << " ms"
              << " | avg:\t " << timing.averageTime() / 1e3 << " us\n";
  }
  std::cout << "============================\n";
}

// ---------------- Scope ----------------
Scope::Scope(const std::string& name)
  : name_(name)
  , start_(std::chrono::high_resolution_clock::now())
{
}

Scope::~Scope()
{
  auto end = std::chrono::high_resolution_clock::now();
  double dt = std::chrono::duration<double, std::nano>(end - start_).count();
  Profiler::instance().addTiming(name_, dt);
}
