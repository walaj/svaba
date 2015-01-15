//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SimplePool - Zero-overhead templated memory pool 
// The design is based on the premise that the allocation
// of the objects have a lifetime the length of the program's
// execution so they are never freed. Should only be used when
// the number of objects is bounded so the pool does not constantly grow.
//
// Not thread-safe.
// 
#ifndef SIMPLEPOOL_H
#define SIMPLEPOOL_H

template<class T>
class SimplePool
{
    public:

        SimplePool()
        {
            size_t bytes_per_object = sizeof(T);
            size_t total_bytes = NUM_OBJECTS * bytes_per_object;
            m_pPool = malloc(total_bytes);
            if(m_pPool == NULL)
            {
                std::cerr << "SimpleStorage failed to allocate " << total_bytes << 
                " bytes for memory pool, exiting\n";
                abort();
            }
            m_capacity = total_bytes;
            m_used = 0;
        }

        ~SimplePool()
        {
            free(m_pPool);
        }
    
        // Return a pointer to the next unused block of memory
        void* alloc()
        {
            assert(m_used < m_capacity);
            void* pNext = (char*)m_pPool + m_used;
            m_used += sizeof(T);
            return pNext;
        }

        // Do not track deallocations
        void dealloc(void* /*ptr*/)
        {
            // does nothing
        }

        bool isFull()
        {
            return m_used >= m_capacity;
        }

    private:

        void* m_pPool;
        size_t m_capacity;
        size_t m_used;
        static const size_t NUM_OBJECTS = 50*1024;
};

#endif
