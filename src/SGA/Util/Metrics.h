///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Metrics - Data structures to record data metrics
// error rates, error distributions, etc
//
#ifndef METRICS_H
#define METRICS_H

#include <map>

struct ErrorCount
{
    int64_t num_samples;
    int64_t num_errors;
};

template<class Key>
class ErrorCountMap
{
    typedef std::map<Key, ErrorCount> DataMap;

    public:
        ErrorCountMap() {}

        //
        void incrementSample(const Key& key)
        {
            ++m_data[key].num_samples;
        }

        //
        void incrementError(const Key& key)
        {
            ++m_data[key].num_errors;
        }

        //
        void write(std::ostream* pWriter, const std::string& leader, const std::string& header)
        {
            *pWriter << leader;
            *pWriter << header << "\tsamples\terrors\tfraction\n";
            for(typename DataMap::iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
            {
                *pWriter << iter->first << "\t" << 
                             iter->second.num_samples << "\t" << 
                             iter->second.num_errors << "\t" << 
                             (double)iter->second.num_errors / iter->second.num_samples << "\n";
            }
        }

    private:

        DataMap m_data;
};

#endif
