/*! \file cmd_line.h
*  \brief simple command line parsing
* \author Hanno Hildenbrandt
*/

#ifndef JC_CMD_LINE_H_INCLUDED
#define JC_CMD_LINE_H_INCLUDED

#include <stdexcept>
#include <cassert>
#include <string.h>
#include <sstream>
#include <utility>
#include <string>


namespace cmd {


  class parse_error : public std::invalid_argument
  {
  public:
    parse_error(const char* msg) : std::invalid_argument(msg)
    {}
  };


  // split argument at '='
  std::pair<std::string, std::string> split_arg(const char* carg)
  {
    const char* s = strchr(carg, '=');
    if (nullptr == s)
    {
      return{ "", "" };
    }
    return{{carg, s}, {s + 1}};
  }


  template <typename T>
  void convert_arg(std::pair<std::string, std::string> const& arg, T& x)
  {
    std::istringstream iss(arg.second);
    if (!(iss >> x))
    {
      throw parse_error((std::string("invalid value for argument ") + arg.first).c_str());
    }
  }


  void parse_cmd_flag(const char* name, bool& val, int argc, const char** argv)
  {
    for (int i = 0; i < argc; ++i)
    {
      if (0 == strcmp(argv[i], name)) { val = true; return; }
    }
  }


  template <typename T, typename P = T>
  void parse_optional_arg(const char* name, T& val, int argc, const char** argv)
  {
    for (int i = 0; i < argc; ++i)
    {
      auto arg = split_arg(argv[i]);
      P pval;
      if (arg.first == name) { convert_arg(arg, pval); val = static_cast<T>(pval);  return; }
    }
  }


  template <typename T, typename P = T>
  void parse_required_arg(const char* name, T& val, int argc, const char** argv)
  {
    int i = 0;
    for (; i < argc; ++i)
    {
      auto arg = split_arg(argv[i]);
      P pval;
      if (arg.first == name) { convert_arg(arg, pval); val = static_cast<T>(pval);  return; }
    }
    throw parse_error(((std::string("missing argument '") + name) + '\'').c_str());
  }

}

#endif
