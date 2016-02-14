#ifndef SNOWMAN_PON_WALKER_H__
#define SNOWMAN_PON_WALKER_H__

#include <pthread.h>
#include <fstream>

#include "SnowTools/BamWalker.h"

#include "DiscordantCluster.h"

class DiscRead: public SnowTools::GenomicRegion
{

 public:
  DiscRead() {}
  DiscRead (int32_t c, uint32_t p1, uint32_t p2, bool s) : SnowTools::GenomicRegion(c, p1, p2, s) {}
    DiscRead (int32_t c, uint32_t p1, bool s, int32_t cm, uint32_t p1m, bool sm, uint16_t mid) : SnowTools::GenomicRegion(c, p1, p1, s) 
    {
      mate = SnowTools::GenomicRegion(cm, p1m, p1m, sm);
      id = mid;
    }

    
  SnowTools::GenomicRegion mate;
  uint16_t id = 0;

  std::string toString() const; 

  std::string toTabString() const; 

  void save(std::ofstream& of);

  void load(std::ifstream& inf) ;


};

typedef std::vector<DiscRead> DiscReadVector;

class PonWalker: public SnowTools::BamWalker {
  
 public:
  void walkBam(std::ofstream& of, pthread_mutex_t* lock);

  PonWalker(const std::string& in, uint16_t mid, SnowTools::BamWalker *o) : SnowTools::BamWalker(in) {
    id = mid;
    m_bw = o;
  }
  PonWalker() : SnowTools::BamWalker() {}


  SnowTools::BamWalker * m_bw;  
 private:
  
  uint16_t id = 0;
  SnowTools::DiscordantClusterMap m_dmap;
  DiscReadVector m_disc;



};

#endif
