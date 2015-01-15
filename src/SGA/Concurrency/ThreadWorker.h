///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ThreadWorker - Generic thread class
// to perform a batch of work. It gets a vector of
// Input items from the master thread, uses Processor to
// perform some operation on the data which returns a value 
// of type Output. A vector of output is swapped back to the master
// thread
//
#ifndef THREADWORKER_H
#define THREADWORKER_H

#include <semaphore.h>
#include "Util.h"

template<class Input, class Output, class Processor>
class ThreadWorker
{
    typedef std::vector<Input> InputVector;
    typedef std::vector<Output> OutputVector;

    public:
        ThreadWorker(sem_t* pReadySem, Processor* pProcessor, const size_t max_items);
        ~ThreadWorker();

        // Exchange the contents of the shared input/output vectors with pInput/pOutput
        void swapBuffers(InputVector& otherInputVector, OutputVector& otherOutputVector);

        // External control functions
        void start();
        void stop();
        bool isReady();

    private:

        // Main work loop
        void run();

        // Main overlap function
        void process(Input& item);
    
        // Thread entry point
        static void* startThread(void* obj);
        
        // Handles
        pthread_t m_thread;

        // External semaphore to post to
        // when the thread is ready to receive data
        sem_t* m_pReadySem;

        // Semaphore the external caller posts to when
        // data is ready to consume
        sem_t m_producedSem;
        
        // Shared data
        pthread_mutex_t m_mutex;
        InputVector m_sharedInputVector;
        OutputVector m_sharedOutputVector;
        Processor* m_pProcessor;

        volatile bool m_stopRequested;
        bool m_isReady;
};

// Implementation
template<class Input, class Output, class Processor>
ThreadWorker<Input, Output, Processor>::ThreadWorker(sem_t* pReadySem, 
                                                     Processor* pProcessor,
                                                     const size_t max_items) :
                                                      m_pReadySem(pReadySem),
                                                      m_pProcessor(pProcessor),
                                                      m_stopRequested(false), 
                                                      m_isReady(false)
{
    m_sharedInputVector.reserve(max_items);
    m_sharedOutputVector.reserve(max_items);

    // Set up semaphores and mutexes
    int ret = sem_init( &m_producedSem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    ret = pthread_mutex_init(&m_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
template<class Input, class Output, class Processor>
ThreadWorker<Input, Output, Processor>::~ThreadWorker()
{
    sem_destroy(&m_producedSem);
    int ret = pthread_mutex_destroy(&m_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}

// Externally-called function from the main thread
// It locks the mutex and swaps the vectors, providing
// data for the work thread to operate on
template<class Input, class Output, class Processor>
void ThreadWorker<Input, Output, Processor>::swapBuffers(InputVector& otherInputVector, OutputVector& otherOutputVector)
{
    pthread_mutex_lock(&m_mutex);
    
    m_isReady = false;
    m_sharedInputVector.swap(otherInputVector);
    m_sharedOutputVector.swap(otherOutputVector);
    pthread_mutex_unlock(&m_mutex);
    sem_post(&m_producedSem);
}

// Externally-called function to start the worker
template<class Input, class Output, class Processor>
void ThreadWorker<Input, Output, Processor>::start()
{
    int ret = pthread_create(&m_thread, 0, &ThreadWorker<Input, Output, Processor>::startThread, this);
    if(ret != 0)
    {
        std::cerr << "Thread creation failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}

// Externally-called function to tell the worker to stop
template<class Input, class Output, class Processor>
void ThreadWorker<Input, Output, Processor>::stop()
{
    // Wait for the work thread to finish
    // with the data its processing then set the stop
    // flag and unblock the work process
    sem_wait(m_pReadySem);
    m_stopRequested = true;
    sem_post(&m_producedSem);
    int ret = pthread_join(m_thread, NULL);
    if(ret != 0)
    {
        std::cerr << "Thread join failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }

}

// 
template<class Input, class Output, class Processor>
bool ThreadWorker<Input, Output, Processor>::isReady()
{
    return m_isReady;
}

// Main worker loop
template<class Input, class Output, class Processor>
void ThreadWorker<Input, Output, Processor>::run()
{
    // Indicate that the thread is ready to receive data
    pthread_mutex_lock(&m_mutex);
    m_isReady = true;
    sem_post(m_pReadySem);
    pthread_mutex_unlock(&m_mutex);

    while(1)
    {
        // Block until there is some data to use
        sem_wait(&m_producedSem);
        
        // If the stop flag is now set, finish
        if(m_stopRequested)
            break;

        // Lock the shared buffer and process the input work
        pthread_mutex_lock(&m_mutex);
        
        assert(m_sharedOutputVector.empty());
        size_t num_items = m_sharedInputVector.size();
        for(size_t i = 0; i < num_items; ++i)
        {
            Output result = m_pProcessor->process(m_sharedInputVector[i]);
            m_sharedOutputVector.push_back(result);    
        }
        
        m_isReady = true;
        pthread_mutex_unlock(&m_mutex);

        // Post to the semaphore
        sem_post(m_pReadySem);
    }
}

// Thread entry point
template<class Input, class Output, class Processor>
void* ThreadWorker<Input, Output, Processor>::startThread(void* obj)
{
    reinterpret_cast<ThreadWorker*>(obj)->run();
    return NULL;
}

#endif
