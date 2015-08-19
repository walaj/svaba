#ifndef GENPON_H
#define GENPON_H

#include "PonWalker.h"
#include <string>
#include <fstream>
#include <pthread.h>
#include "workqueue.h"
#include "SnowTools/GenomicRegionCollection.h"

void parsePONOptions(int argc, char** argv);
void runGeneratePONfromVCF(int argc, char** argv);
void runGeneratePON(int argc, char** argv);
void combinePON();
void verifyPON();


class PONThread {

  public:
  PONThread() : m_tid(0), m_running(0), m_detached(0) {}
  
  ~PONThread()
  {
    if (m_running == 1 && m_detached == 0) {
      pthread_detach(m_tid);
    }
    if (m_running == 1) {
      pthread_cancel(m_tid);
    }
  }

  static void* runThread(void* arg)
  {
    return ((PONThread*)arg)->run();
  }
 
  int start()
  {
    int result = pthread_create(&m_tid, NULL, runThread, this);
    if (result == 0) {
      m_running = 1;
    }
    return result;
  }

  int join()
  {
    int result = -1;
    if (m_running == 1) {
      result = pthread_join(m_tid, NULL);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }
 
  int detach()
  {
    int result = -1;
    if (m_running == 1 && m_detached == 0) {
      result = pthread_detach(m_tid);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }

  pthread_t self() {
    return m_tid;
  }

  virtual void* run() = 0;
 
private:
  pthread_t  m_tid;
  int        m_running;
  int        m_detached;
};

template <class T>
  class PONConsumerThread : public PONThread {
 
public:

  PONConsumerThread(wqueue<T*>& queue, bool verbose) : m_queue(queue), m_verbose(verbose) {}
 
  void* run() {
    // Remove 1 item at a time and process it. Blocks if no items are 
    // available to process.
    for (int i = 0;; i++) {
      //if (m_verbose)
	//printf("thread %lu, loop %d - waiting for item...\n", 
	//     (long unsigned int)self(), i);
      T* item = (T*)m_queue.remove();
      item->run();
      //m_output->push_back(item->output());
      delete item;
      if (m_queue.size() == 0)
        return NULL;
    }
    return NULL;
  }

  //std::vector<O*> getThreadOutput() {return m_output;}

 private: 
  wqueue<T*>& m_queue;
  bool m_verbose;

};

class PONWorkItem {

private:

  std::string m_bam;
  std::ofstream * m_of;
  pthread_mutex_t * m_lock;
  size_t m_id;
  SnowTools::MiniRulesCollection * m_mr;
  SnowTools::GRC * m_fr;
  SnowTools::BamWalker *m_out;
  
public:

  PONWorkItem(const std::string& bam, std::ofstream * of, pthread_mutex_t * lock, size_t id, 
	      SnowTools::MiniRulesCollection* mr, SnowTools::GRC *fr, SnowTools::BamWalker *bw)
    : m_bam(bam), m_of(of), m_lock(lock), m_id(id), m_mr(mr), m_fr(fr), m_out(bw)
  {}

  /// Destroy this work item
  ~PONWorkItem() { }


  /** Kick off all of the swapping and annealing
   * @return true if allPONs completed
   */
  bool run() 
  { 
    PonWalker pw(m_bam, m_id, m_out);
    if (m_fr->size()) 
      pw.setBamWalkerRegions(m_fr->asGenomicRegionVector());
    pthread_mutex_lock(m_lock);  
    std::cerr << "...walking BAM to generate Panel of Normal data for BAM: "  << m_bam << std::endl;
    pthread_mutex_unlock(m_lock);  
    pw.SetMiniRulesCollection(*m_mr);
    pw.walkBam(*m_of, m_lock);
    return true;
     
  }

};


#endif

/*  LocalWords:  arg
 */
