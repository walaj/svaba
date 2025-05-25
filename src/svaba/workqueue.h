#pragma once

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <vector>
#include <memory>
#include <atomic>
#include <map>
#include "svabaThreadUnit.h"
#include "SeqLib/RefGenome.h"

// A simple thread-safe queue for unique_ptr jobs
template<typename Job>
class WorkQueue {
public:
  void push(std::unique_ptr<Job> job) {
    {
      std::lock_guard<std::mutex> lk(mtx_);
      queue_.push_back(std::move(job));
    }
    cv_.notify_one();
  }

  std::unique_ptr<Job> pop() {
    std::unique_lock<std::mutex> lk(mtx_);
    cv_.wait(lk, [this]{ return !queue_.empty(); });
    auto job = std::move(queue_.front());
    queue_.pop_front();
    return job;
  }

  bool empty() const {
    std::lock_guard<std::mutex> lk(mtx_);
    return queue_.empty();
  }

private:
  mutable std::mutex mtx_;
  std::condition_variable cv_;
  std::deque<std::unique_ptr<Job>> queue_;
};

// A pool of worker threads that pull WorkItems from the queue
template<typename WorkItem>
class ThreadPool {
public:
  using JobPtr = std::unique_ptr<WorkItem>;

  ThreadPool(size_t numThreads,
             const std::string& refIndex,
             const std::string& virIndex,
             const std::map<std::string,std::string>& bamFiles)
    : running_(true)
  {
    for (size_t i = 0; i < numThreads; ++i) {
      workers_.emplace_back([=]{
        // Thread-local setup:
        svabaThreadUnit wu;
        wu.ref_genome = new SeqLib::RefGenome();
        wu.ref_genome->LoadIndex(refIndex);
        if (!virIndex.empty()) {
          wu.vir_genome = new SeqLib::RefGenome();
          wu.vir_genome->LoadIndex(virIndex);
        }
        for (auto& p : bamFiles) {
          wu.walkers[p.first].Open(p.second);
          wu.walkers[p.first].prefix = p.first;
        }

        // Main loop:
        while (running_) {
          auto job = queue_.pop();
          if (!job) break;            // nullptr sentinel to exit
          job->run(wu, std::this_thread::get_id());
        }
      });
    }
  }

  // Submit a new WorkItem to the queue
  void enqueue(JobPtr job) {
    queue_.push(std::move(job));
  }

  // Gracefully shut down the pool: send sentinel and join
  void shutdown() {
    running_ = false;
    // push one null job per thread to unblock
    for (size_t i = 0; i < workers_.size(); ++i)
      queue_.push(nullptr);
    for (auto& t : workers_)
      t.join();
  }

private:
  std::atomic<bool> running_;
  WorkQueue<WorkItem> queue_;
  std::vector<std::thread> workers_;
};

/*
#ifndef WORKQUEUE_SVABA_H
#define WORKQUEUE_SVABA_H

#include <pthread.h>
#include <list>

#include "svabaThreadUnit.h"
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

class svabaThread {

  public:
  svabaThread() : m_tid(0), m_running(0), m_detached(0) {}
  
  ~svabaThread()
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
    return ((svabaThread*)arg)->run();
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
  class ConsumerThread : public svabaThread {
 
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

  svabaThreadUnit wu;

 private: 
  wqueue<T*>& m_queue;
  bool m_verbose;

};

#endif
*/
