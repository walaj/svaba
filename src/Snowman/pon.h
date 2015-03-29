#ifndef SNOWMAN_PON_GEN_H
#define SNOWMAN_PON_GEN_H

#include <cstdlib>
#include <string>
#include <unordered_map>


typedef std::unordered_map<std::string, size_t> PON;

void runPON(int argc, char** argv);
void parsePONOptions(int argc, char** argv);

#endif


