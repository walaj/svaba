#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include "GenomicRegion.h"

bool runExtractor(int argc, char** argv);
void parseExtractOptions(int argc, char** argv);
bool parseRegionFile(GenomicRegionVector &gr);

#endif
