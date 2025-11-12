#pragma once
#include <format>
#include <string>
#include <utility>

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
#else
constexpr void Init(LoggerType type) { }
constexpr void Close() { }
constexpr void Debug(const std::string& string) { }
constexpr void Warning(const std::string& string) { }
constexpr void Error(const std::string& string) { }
#endif
#define DebugF(str, ...) LOG::Debug(std::format(str, __VA_ARGS__))
#define WarningF(str, ...) LOG::Warning(std::format(str, __VA_ARGS__))
#define ErrorF(str, ...) LOG::Error(std::format(str, __VA_ARGS__))
}
