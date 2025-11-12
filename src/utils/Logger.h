#pragma once
#ifdef DEBUG
#include <format>
#include <utility>
#endif
#include <string>

namespace LOG {
enum class LoggerType
{
  STDOUT,
  STDERR,
  FILE
};
#ifdef DEBUG
void Init(LoggerType type);
void Close();
void Debug(const std::string& string);
void Warning(const std::string& string);
void Error(const std::string& string);
#define DebugF(str, ...) LOG::Debug(std::format(str, __VA_ARGS__))
#define WarningF(str, ...) LOG::Warning(std::format(str, __VA_ARGS__))
#define ErrorF(str, ...) LOG::Error(std::format(str, __VA_ARGS__))
#else
inline void Init(LoggerType type) { }
inline void Close() { }
inline void Debug(const std::string& string) { }
inline void Warning(const std::string& string) { }
inline void Error(const std::string& string) { }
#define DebugF(str, ...) 0
#define WarningF(str, ...) 0
#define ErrorF(str, ...) 0
#endif
}
