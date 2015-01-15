#include <pthread.h>

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

  ConsumerThread(wqueue<T*>& queue, bool verbose) : m_queue(queue), m_verbose(verbose) {}
 
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
