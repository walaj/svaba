#ifndef SNOWTOOLS_UTIL_H
#define SNOWTOOLS_UTIL_H

#include <time.h>
#include <string>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace BamTools;
using namespace std;

void displayRuntime(const timespec start);
string getSamHeader(string bamfile, SamHeader &sam);
void getRefVector(string bamfile, RefVector &ref);

#endif
