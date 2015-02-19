#ifndef SNOW_CLUSTER_H
#define SNOW_CLUSTER_H

#include "SVBamReader.h"
#include "GenomicRegion.h"
#include <sstream>
#include <iostream>
#include <unordered_map>

typedef unordered_map<string, string> RMap;
typedef unordered_map<string, size_t> GMap;

namespace ClusterReads {
  void clusterReads(BamAlignmentVector &bav, GenomicRegionVector &grv, 
		    RMap &rmap, string orientation, int isize, 
		    int cluster_buffer);
  void finalizeCluster(GenomicRegionVector &grv, RMap &rmap, int pos, GenomicRegion anc, GenomicRegion par, int dtcount, int dncount);
  bool validDiscordant(const BamAlignment &it, int isize, string orientation);
}

#endif
