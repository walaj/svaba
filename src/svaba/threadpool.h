// threadpool.h
#pragma once

#include <deque>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <vector>
#include <atomic>

#include "svabaOutputWriter.h"

// A simple thread safe queue for unique_ptr jobs
template<typename Job>
class WorkQueue {
public:
  void push(std::unique_ptr<Job> job) {
    {
      std::lock_guard lk(mtx_);
      queue_.push_back(std::move(job));
    }
    cv_.notify_one();
  }

  // pop one job; blocks until one is available or a nullptr sentinel is received
  std::unique_ptr<Job> pop() {
    std::unique_lock lk(mtx_);
    cv_.wait(lk, [&]{ return !queue_.empty(); });
    auto job = std::move(queue_.front());
    queue_.pop_front();
    return job;
  }

  // send one nullptr for each worker so they can exit
  void shutdown(size_t numThreads) {
    for(size_t i=0;i<numThreads;++i)
      push(nullptr);
  }

private:
  std::mutex                mtx_;
  std::condition_variable   cv_;
  std::deque<std::unique_ptr<Job>> queue_;
};

// A fixed-size thread pool that takes Job functors
template<typename WorkItem>
class ThreadPool {
public:
  ThreadPool(size_t nThreads,
             const std::string& ref,
             const std::map<std::string,std::string>& bamFiles,
	     SvabaOutputWriter& output_writer)
    : queue_(), workers_(), flushMutex_(), output_writer_(output_writer)
  {
    workers_.reserve(nThreads);
    for(size_t i=0;i<nThreads;++i){
      workers_.emplace_back([=,&bamFiles,&flushMutex=this->flushMutex_](){
	  
        // per-thread setup - each thread gets its own FASTA read and BAM readers
	// these should not be shared across threads, even if using const functions only
        svabaThreadUnit unit;
        unit.ref_genome = std::make_unique<SeqLib::RefGenome>();
        unit.ref_genome->LoadIndex(ref);
	
	// open the BAM files
        for(auto &p : bamFiles){
          unit.walkers[p.first].Open(p.second);
          unit.walkers[p.first].prefix = p.first;
        }

        // consume jobs
        while(auto job = queue_.pop()){
          (*job)(unit, std::hash<std::thread::id>{}(std::this_thread::get_id()));
        }

        // final flush
	output_writer_.writeUnit(unit); // this does the flush and mutex in it
	});
    }
  }

  // submit a new work item
  void submit(std::unique_ptr<WorkItem> job){
    queue_.push(std::move(job));
  }

  // tell threads no more work, then join
  void shutdown(){
    queue_.shutdown(workers_.size());
    for(auto &t : workers_) t.join();
  }

private:
  WorkQueue<WorkItem>      queue_;
  std::vector<std::thread> workers_;
  std::mutex               flushMutex_;
  SvabaOutputWriter&       output_writer_; 

};
