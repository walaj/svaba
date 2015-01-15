///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Implementation of a multikey quicksort worker thread
//
#include <pthread.h>
#include <semaphore.h>
#include "mkqs.h"

// 
template<typename T>
struct MkqsJob
{
    MkqsJob(T* p, int num, int d) : pData(p), n(num), depth(d) {}
    T* pData;
    int n;
    int depth;
};

//
template<typename T, class PrimarySorter, class FinalSorter>
class MkqsThread
{
    typedef MkqsJob<T> Job;
    typedef std::queue<Job> JobQueue;

    public:
        MkqsThread(int id, JobQueue* pQueue, pthread_mutex_t* pQueueMutex, 
                   sem_t* pQueueSem, sem_t* pDoneSem, int thresholdSize,
                   const PrimarySorter* pPrimarySorter, 
                   const FinalSorter* pFinalSorter) : m_id(id), 
                                                      m_pQueue(pQueue), 
                                                      m_pQueueMutex(pQueueMutex), 
                                                      m_pQueueSem(pQueueSem),
                                                      m_pDoneSem(pDoneSem),
                                                      m_thresholdSize(thresholdSize),
                                                      m_pPrimary(pPrimarySorter), 
                                                      m_pFinal(pFinalSorter),
                                                      m_stopRequested(false),
                                                      m_numProcessed(0) {}
        ~MkqsThread();

        void start();
        void stop();
        void join();

        static void* startThread(void* obj);

    private:

        void run();
        void process(Job& job);

        // Data
        int m_id;
        JobQueue* m_pQueue; // shared
        pthread_mutex_t* m_pQueueMutex; // shared
        sem_t* m_pQueueSem; // shared
        sem_t* m_pDoneSem; // shared
        
        int m_thresholdSize;
        const PrimarySorter* m_pPrimary;
        const FinalSorter* m_pFinal;

        pthread_t m_thread;
        volatile bool m_stopRequested;
        int m_numProcessed;
};

//
template<typename T, class PrimarySorter, class FinalSorter>
MkqsThread<T, PrimarySorter, FinalSorter>::~MkqsThread()
{
}

// Start thread
template<typename T, class PrimarySorter, class FinalSorter>
void MkqsThread<T, PrimarySorter, FinalSorter>::start()
{
    int ret = pthread_create(&m_thread, 0, &MkqsThread<T, PrimarySorter, FinalSorter>::startThread, this);
    if(ret != 0)
    {
        std::cerr << "Thread creation failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Stop thread
template<typename T, class PrimarySorter, class FinalSorter>
void MkqsThread<T, PrimarySorter, FinalSorter>::stop()
{
    m_stopRequested = true;
}

// Called from the external main function, joins the thread to the main on exit
template<typename T, class PrimarySorter, class FinalSorter>
void MkqsThread<T, PrimarySorter, FinalSorter>::join()
{
    int ret = pthread_join(m_thread, NULL);
    if(ret != 0)
    {
        std::cerr << "Thread join failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Run thread
template<typename T, class PrimarySorter, class FinalSorter>
void MkqsThread<T, PrimarySorter, FinalSorter>::run()
{
    while(1)
    {
        sem_post(m_pDoneSem);
        sem_wait(m_pQueueSem);
        
        // Exit if the thread was stopped
        if(m_stopRequested)
        {
            pthread_exit(NULL);
        }

        // Take an item from the queue and process it
        pthread_mutex_lock(m_pQueueMutex);
        assert(!m_pQueue->empty());

        // Decrement the done count. Since this is done
        // while the queue is not empty and locked, the master
        // thread can never see an empty queue before the
        // done semaphore is decremented. This ensures all threads
        // will work to completion
        sem_wait(m_pDoneSem);

        Job job = m_pQueue->front();
        m_pQueue->pop();
        pthread_mutex_unlock(m_pQueueMutex);

        // Process the item using either the parallel algorithm (which subdivides the job further)
        // or the serial algorithm (which doesn't subdivide)
        if(job.n > m_thresholdSize)
        {
            parallel_mkqs_process(job, m_pQueue, m_pQueueMutex, m_pQueueSem, *m_pPrimary, *m_pFinal);
        }
        else
        {
            mkqs2(job.pData, job.n, job.depth, *m_pPrimary, *m_pFinal);
        }
        m_numProcessed += 1;
    }
}

// Thread entry point
template<typename T, class PrimarySorter, class FinalSorter>
void* MkqsThread<T, PrimarySorter, FinalSorter>::startThread(void* obj)
{
    reinterpret_cast<MkqsThread*>(obj)->run();
    return NULL;
}
