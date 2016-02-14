#include "PonWalker.h"

using SnowTools::BamRead;
using SnowTools::BamReadVector;

std::string DiscRead::toString() const {
  
  std::stringstream out;
  out << "ID " << id << " " << (chr+1) << ":" << SnowTools::AddCommas<int>(pos1) << "-->" << 
    (mate.chr+1) << ":" << SnowTools::AddCommas<int>(mate.pos1);
  return out.str();
}

std::string DiscRead::toTabString() const {
  
  std::stringstream out;
  out << id << "\t" << (chr+1) << "\t" << pos1 << "\t" << 
    (mate.chr+1) << "\t" << mate.pos1;
  return out.str();
}


void DiscRead::save(std::ofstream& of) {

  char c1 = chr;
  char cm = mate.chr;
  of.write((char*)&id, sizeof(uint16_t)); //sizeof(c1)); 
  of.write(&c1, sizeof(char)); //sizeof(c1)); 
  of.write((char*)&pos1, sizeof(int32_t)); 
  of.write(&cm, sizeof(uint8_t)); 
  of.write((char*)&(mate.pos1), sizeof(int32_t)); 

}

void DiscRead::load(std::ifstream& inf) 
{ 
  
  char c1, cm;
  inf.read((char*)&id, sizeof(uint16_t));   
  inf.read(&c1, sizeof(char)); 
  inf.read((char*)&pos1, sizeof(int32_t)); 
  inf.read(&cm, sizeof(char)); 
  inf.read((char*)&mate.pos1, sizeof(int32_t)); 

  chr = c1;
  mate.chr = cm;
  pos2 = pos1;
  mate.pos2 = mate.pos1;

} 
  
    


void PonWalker::walkBam(std::ofstream& of, pthread_mutex_t* lock) {

  size_t CHUNKSIZE = 1000000;

  BamRead r;

  bool rule_pass;
  
  SnowTools::BamReadVector buff;
  SnowTools::GenomicRegion region;
  size_t countr = 0;

  std::string sid = std::to_string(id);
  while(GetNextRead(r, rule_pass)) {

    ++countr;

    if (countr % 1000000 == 0)
      std::cerr << "...read " << SnowTools::AddCommas<size_t>(countr) <<  " at " << (r.ChrID() + 1) << ":" << SnowTools::AddCommas<int>(r.Position()) << std::endl;

    if (rule_pass) {
      r.clearSeqQualAndTags();
      r.SetQname(sid);
      buff.push_back(r);
    }
    
    //if (rule_pass)
    //  m_disc.push_back(DiscRead(r.ChrID(), r.Position(), !r.ReverseFlag(), r.MateChrID(), r.MatePosition(), !r.MateReverseFlag(), id));
    
    if (m_disc.size() > CHUNKSIZE) {
      pthread_mutex_lock(lock);  
      std::cerr << "...writing to PON BAM for BAM id: " << id << std::endl; 
      for (auto& i : buff)
	m_bw->writeAlignment(i);
      pthread_mutex_unlock(lock);  
      buff.clear();
    }
  }

  // finish it
  pthread_mutex_lock(lock);  
  std::cerr << "...writing final for BAM " << id << std::endl; 
  for (auto& i : buff)
    m_bw->writeAlignment(i);
    //i.save(of);
  pthread_mutex_unlock(lock);  
  m_disc.clear();
  
  
}
