#ifndef __CL_ARG_PARSING_H__
#define __CL_ARG_PARSING_H__

#include <string>
#include <sstream>
#include <iostream>

//========================================================
// code to read input arguments
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}
// Special version for string without the conversion
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}

#endif
