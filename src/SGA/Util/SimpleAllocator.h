//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SimpleAllocator - High-level manager of SimplePool
// memory pools. See SimplePool.h for description of allocation
// strategy
//
#ifndef SIMPLEALLOCATOR_H
#define SIMPLEALLOCATOR_H

#include <list>
#include "SimplePool.h"

template<class T>
class SimpleAllocator
{
    typedef SimplePool<T> StorageType;
    typedef std::list<StorageType* > StorageList;

    public:
        SimpleAllocator() {}

        ~SimpleAllocator()
        {
            for(typename StorageList::iterator iter = m_pPoolList.begin(); iter != m_pPoolList.end(); ++iter)
            {
                delete *iter;
            }
            m_pPoolList.clear();
        }

        void* alloc()
        {
            if(m_pPoolList.empty() || m_pPoolList.back()->isFull())
            {
                // new storage must be allocated
                m_pPoolList.push_back(new StorageType);
            }

            // allocate from the last pool
            return m_pPoolList.back()->alloc();
        }

        void dealloc(void* /*ptr*/)
        {
            // deallocation not tracked in this strategy
        }

    private:

        StorageList m_pPoolList;
};

#endif
