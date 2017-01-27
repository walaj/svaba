#ifndef WORKQUEUE_SNOW_H
#define WORKQUEUE_SNOW_H

#include <pthread.h>
#include <list>

#include "svabaWorkUnit.h"
#include "SeqLib/RefGenome.h"

typedef std::map<std::string, svabaBamWalker> WalkerMap;

template <typename T> class wqueue
{ 

  public:
  wqueue() {
    pthread_mutex_init(&m_mutex, NULL);
    pthread_cond_init(&m_condv, NULL);
  }

  ~wqueue() {
    pthread_mutex_destroy(&m_mutex);
    pthread_cond_destroy(&m_condv);
  }

  void add(T item) {
    pthread_mutex_lock(&m_mutex);
    m_queue.push_back(item);
    pthread_cond_signal(&m_condv);
    pthread_mutex_unlock(&m_mutex);
  }

  T remove() {
    pthread_mutex_lock(&m_mutex);
    while (m_queue.size() == 0) {
      pthread_cond_wait(&m_condv, &m_mutex);
    }
    T item = m_queue.front();
    m_queue.pop_front();
    pthread_mutex_unlock(&m_mutex);
    return item;
  }

  int size() {
    pthread_mutex_lock(&m_mutex);
    int size = m_queue.size();
    pthread_mutex_unlock(&m_mutex);
    return size;
  }

  std::list<T>   m_queue;
  pthread_mutex_t m_mutex;
  pthread_cond_t  m_condv;

};

class SnowThread {

  public:
  SnowThread() : m_tid(0), m_running(0), m_detached(0) {}
  
  ~SnowThread()
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
    return ((SnowThread*)arg)->run();
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
  class ConsumerThread : public SnowThread {
 
public:

 ConsumerThread(wqueue<T*>& queue, bool verbose, 
		const std::string& ref, const std::string& vir,
		const std::map<std::string, std::string>& bams) : m_queue(queue), m_verbose(verbose) {

    // load the reference genomce
    if (m_verbose)
      std::cerr << "\tOpening ref genome for thread " << self() << std::endl;
    wu.ref_genome = new SeqLib::RefGenome();
    wu.ref_genome->LoadIndex(ref);

    // load the viral genome
    if (!vir.empty()) {
      if (m_verbose)
	std::cerr << "\tOpening vir genome for thread " << self() << std::endl;
      wu.vir_genome = new SeqLib::RefGenome();
      wu.vir_genome->LoadIndex(vir);
    } 

    // open the bams for this thread
    if (m_verbose)
      std::cerr << "\tOpening BAMs for thread " << self() << std::endl;
    for (auto& b : bams) {
      wu.walkers[b.first] = svabaBamWalker();
      wu.walkers[b.first].Open(b.second);
      wu.walkers[b.first].prefix = b.first;
    }
    

  }
 
  void* run() {
    // Remove 1 item at a time and process it. Blocks if no items are 
    // available to process.
    for (int i = 0;; i++) {
      //if (m_verbose)
	//printf("thread %lu, loop %d - waiting for item...\n", 
	//     (long unsigned int)self(), i);
      T* item = (T*)m_queue.remove();
      item->run(wu, (long unsigned)self()); 
      delete item;
      if (m_queue.size() == 0)
        return NULL;
    }
    return NULL;
  }

  svabaWorkUnit wu;

 private: 
  wqueue<T*>& m_queue;
  bool m_verbose;

};

#endif
