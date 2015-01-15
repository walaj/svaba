//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Timer - Simple object to that prints the wallclock 
// length of its lifetime
//
#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <sys/time.h>

class Timer
{
    public:
        Timer(std::string s, bool silent = false) : m_desc(s), m_silent(silent)
        {
            reset();
        }

        ~Timer() 
        {
            if(!m_silent) 
                fprintf(stderr, "[timer - %s] wall clock: %.2lfs CPU: %.2lfs\n", m_desc.c_str(), getElapsedWallTime(), getElapsedCPUTime()); 
        }

        double getElapsedWallTime() const 
        { 
            timeval now;
            gettimeofday(&now, NULL);
            return (now.tv_sec - m_wallStart.tv_sec) + (double(now.tv_usec - m_wallStart.tv_usec) / 1000000);
        }

        double getElapsedCPUTime() const
        {
            double now = clock();
            return (now - m_cpuStart) / CLOCKS_PER_SEC;
        }

        void reset() { gettimeofday(&m_wallStart, NULL); m_cpuStart = clock(); }

    private:
        std::string m_desc;
        
        // Track the wall-clock and CPU times
        // CPU time includes all threads
        timeval m_wallStart;
        double m_cpuStart;

        bool m_silent;
};

#endif
