//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// mkqs - multikey quicksort
//
// Perform a ternary quicksort of strings as described in
// Bentley and Sedgewick, 1997
//
// Example code was downloaded from http://www.cs.princeton.edu/~rs/strings/demo.c
// Modified by JTS to take in a comparator and use a generic type
//

#ifndef MKQS_H
#define MKQS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <queue>
//#include "MkqsThread.h"

#define mkqs_swap(a, b) { T tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

// Swap [i..i+n] and [j..j+n] in x
template<typename T>
inline void vecswap(int i, int j, int n, T* x)
{   
    while (n-- > 0) 
    {
        mkqs_swap(i, j);
        i++;
        j++;
    }
}

#define mkqs_swap2(a, b) { T t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(p) (primarySorter.getChar(*(p), depth))
#define elem2char(e, d) (primarySorter.getChar((e), (d)))

template<typename T>
inline void vecswap2(T* a, T* b, int n)
{   while (n-- > 0) 
    {
        T t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

template<typename T, typename PrimarySorter>
T* med3func(T* a, T* b, T* c, int depth, const PrimarySorter& primarySorter)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, depth, primarySorter)
template<typename T, typename PrimarySorter, typename FinalSorter>
inline void inssort(T* a, int n, int d, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
    T *pi, *pj, s, t;
    for (pi = a + 1; --n > 0; pi++)
    {
        for (pj = pi; pj > a; pj--) 
        {
            // Inline strcmp: break if *(pj-1) <= *pj
            T elem_s = *(pj - 1);
            T elem_t = *pj;
            const char* s = primarySorter.getChrPtr(elem_s);
            const char* t = primarySorter.getChrPtr(elem_t);

            for (s=s+d, t=t+d; *s==*t && *s!=0; s++, t++)
                ;
            if (*s < *t || (*s == *t && finalSorter(elem_s, elem_t)))
                break;
            mkqs_swap2(pj, pj-1);
        }
    }
}

// 
template<typename T, typename PrimarySorter, typename FinalSorter>
void mkqs2(T* a, int n, int depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
    int r, partval;
    T *pa, *pb, *pc, *pd, *pm, *pn, t;
   
    if (n < 10) 
    {
        inssort(a, n, depth, primarySorter, finalSorter);
        return;
    }

    pm = a + (n/2);
    pn = a + (n-1);

    int mid_idx = rand() % n;

    pm = &a[mid_idx];
    mkqs_swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) 
    {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) 
        {
            if (r == 0) { mkqs_swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) 
        {
            if (r == 0) { mkqs_swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        mkqs_swap2(pb, pc);
        pb++;
        pc--;
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
    if ((r = pb-pa) > 1)
        mkqs2(a, r, depth, primarySorter, finalSorter);
    if (ptr2char(a + r) != 0)
        mkqs2(a + r, pa-a + pn-pd-1, depth+1, primarySorter, finalSorter);
    else
    {
        int n2 = pa - a + pn - pd - 1;
        std::sort(a + r, a + r + n2, finalSorter);
    }
    if ((r = pd-pc) > 1)
        mkqs2(a + n-r, r, depth, primarySorter, finalSorter);
}

// Parallel multikey quicksort. It performs mkqs but will
// subdivide the array to sort into sub jobs which can be sorted using threads.
/*template<typename T, typename PrimarySorter, typename FinalSorter>
void parallel_mkqs(T* pData, int n, int numThreads, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{
    typedef MkqsJob<T> Job;
    typedef std::queue<Job> JobQueue;
    Job initialJob(pData, n, 0);
    JobQueue queue;
    queue.push(initialJob);
    
    // Create the mutex that guards the queue
    pthread_mutex_t queue_mutex;
    int ret = pthread_mutex_init(&queue_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Calculate the threshold size for performing serial continuation of the sort. Once the chunks 
    // are below this size, it is better to not subdivide the problem into smaller chunks
    // to avoid the overhead of locking, adding to the queue, etc. 
    int threshold_size = n / numThreads;

    // Create the semaphore used to signal that data is ready to be processed
    // Initial value is 1 as there is one item on the queue to start
    sem_t queue_sem;
    ret = sem_init( &queue_sem, PTHREAD_PROCESS_PRIVATE, 1 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    // Create the semaphore indicating a thread is finished working.
    // This semaphore is incremented by the threads in their run loop
    // before waiting for queue items. If the thread takes an item,
    // this semaphore is decremented. The main thread checks
    // this semaphore after confirming the queue is empty. If this 
    // semaphore value equals the total number of threads, 
    // no more work can remain the the threads are cleaned up
    sem_t done_sem;
    ret = sem_init( &done_sem, PTHREAD_PROCESS_PRIVATE, 0 );
    if(ret != 0)
    {
        std::cerr << "Semaphore initialization failed with error " << ret << "\n";
        std::cerr << "You are probably running on OSX which does not provide unnamed semaphores\n";
        exit(EXIT_FAILURE);
    }

    // Create and start the threads
    MkqsThread<T, PrimarySorter, FinalSorter>* threads[numThreads];
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i] = new MkqsThread<T, PrimarySorter, FinalSorter>(i, &queue, &queue_mutex, &queue_sem, &done_sem, threshold_size, &primarySorter, &finalSorter);   
        threads[i]->start();
    }

    // Check for the end condition
    bool done = false;
    while(!done)
    {
      //        sleep(1);

        // Check if the queue is empty
        // If it is and all threads are finished working (all have posted to
        // the done semaphore), then the threads can be cleaned up.
        pthread_mutex_lock(&queue_mutex);
        if(queue.empty())
        {
            int semval;
            sem_getvalue(&done_sem, &semval);
            if(semval == numThreads)
                done = true;
        }
        pthread_mutex_unlock(&queue_mutex);
    }

    // Signal all the threads to stop, then post to the semaphore they are waiting on
    // All threads will pick up the stop request after the posts and call pthread exit
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i]->stop();
    }

    for(int i = 0; i < numThreads; ++i)
    {
        sem_post(&queue_sem);
    }

    // Join and destroy the threads
    for(int i = 0; i < numThreads; ++i)
    {
        threads[i]->join();
        delete threads[i];
    }

    // Destroy the semaphore
    sem_destroy(&queue_sem);
    sem_destroy(&done_sem);

    // Destroy the queue mutex
    ret = pthread_mutex_destroy(&queue_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}
*/
//
// Perform a partial sort of the data using the mkqs algorithm
// Iterative sort jobs are created and added to pQueue which is
// protected by pQueueMutex. After addition, pQueueSem is updated.
//
/*
template<typename T, class PrimarySorter, class FinalSorter>
void parallel_mkqs_process(MkqsJob<T>& job, 
                           std::queue<MkqsJob<T> >* pQueue, 
                           pthread_mutex_t* pQueueMutex, 
                           sem_t* pQueueSem,
                           const PrimarySorter& primarySorter, 
                           const FinalSorter& finalSorter)
{
    T* a = job.pData;
    int n = job.n;
    int depth = job.depth;
    
    int r, partval;
    T *pa, *pb, *pc, *pd, *pm, *pn, t;
    
    if(n < 10) 
    {
        inssort(a, n, depth, primarySorter, finalSorter);
        return;
    }
    
    pm = a + (n/2);
    pn = a + (n-1);

    int mid_idx = rand() % n;

    pm = &a[mid_idx];
    mkqs_swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) 
    {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) 
        {
            if (r == 0) { mkqs_swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) 
        {
            if (r == 0) { mkqs_swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        mkqs_swap2(pb, pc);
        pb++;
        pc--;
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

    // Lock the queue and push new items if necessary
    // If new items are added to the queue, the semaphore is posted to
    pthread_mutex_lock(pQueueMutex);

    if ((r = pb-pa) > 1)
    {
        MkqsJob<T> job(a, r, depth);
        pQueue->push(job);
        sem_post(pQueueSem);
    }
    
    if (ptr2char(a + r) != 0)
    {
        MkqsJob<T> job(a + r, pa-a + pn-pd-1, depth + 1);
        pQueue->push(job);
        sem_post(pQueueSem);
    }
    else
    {
        // Finalize the sort
        int n2 = pa - a + pn - pd - 1;
        std::sort(a + r, a + r + n2, finalSorter);
    }

    if ((r = pd-pc) > 1)
    {
        MkqsJob<T> job(a + n-r, r, depth);
        pQueue->push(job);
        sem_post(pQueueSem);
    }

    // Unlock the mutex
    pthread_mutex_unlock(pQueueMutex);
}
*/
#endif

