#include "Logger.h"
#include <fstream>
#include <iostream>
#include <ostream>

namespace LOG {
#ifdef DEBUG
static std::ostream* outstream;
static std::ofstream* outfstream;
static LoggerType loggerType;
void Init(LoggerType type)
{
  loggerType = type;
  switch (type)
  {
  case LoggerType::STDOUT:
  default:
    outstream = &std::cout;
    break;
  case LoggerType::STDERR:
    outstream = &std::cerr;
    break;
  case LoggerType::FILE:
    outfstream = new std::ofstream();
    outfstream->open("Log.txt");
    outstream = outfstream;
    break;
  }
}
void Close()
{
  switch (loggerType)
  {
  case LoggerType::STDOUT:
  case LoggerType::STDERR:
  default:
    break;
  case LoggerType::FILE:
    outfstream->close();
    break;
  }
}
void Debug(const std::string& string)
{
  if (outstream)
    *outstream << string << std::endl;
}
void Warning(const std::string& string)
{
  if (outstream)
    *outstream << "WARNING: " << string << std::endl;
}
void Error(const std::string& string)
{
  if (outstream)
    *outstream << "Error: " << string << std::endl;
}
#endif
}
