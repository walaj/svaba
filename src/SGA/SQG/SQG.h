//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SQG - Common definitions and functions for parsing 
// assembly and sequence graph files
//
#ifndef SQG_H
#define SQG_H

#include <vector>
#include <string>
#include "Util.h"

namespace SQG
{
    //
    const char FIELD_SEP = '\t';
    const char TAG_SEP = ':';

    // getTypeCode must be defined for every type that TagValue can take
    static inline char getTypeCode(int)                 { return 'i'; }
    static inline char getTypeCode(char)                { return 'A'; }
    static inline char getTypeCode(const std::string&)  { return 'Z'; }
    static inline char getTypeCode(float)               { return 'f'; }


    // Tokenize functions
    StringVector tokenizeRecord(const std::string& record);
    StringVector tokenizeTagValue(const std::string& tagValue);

    // Two-state TagValue that allows the data value to be not set
    template<typename T>
    class TagValue
    {
        public:
            TagValue() : m_isInitialized(false) {}
            TagValue(const T& v) : m_isInitialized(true), m_value(v) {}

            T get() const
            {
                assert(m_isInitialized);
                return m_value;
            }

            void set(T v)
            {
                m_isInitialized = true;
                m_value = v;
            }

            bool isInitialized() const
            {
                return m_isInitialized;
            }

            std::string toTagString(const char* tag)
            {
                std::stringstream ss;
                char type_code = getTypeCode(m_value);
                ss << tag << TAG_SEP << type_code << TAG_SEP << m_value;
                return ss.str();
            }

            void fromString(const std::string& str)
            {
                StringVector tokens = tokenizeTagValue(str);
                if(tokens.size() != 3)
                {
                    std::cerr << "Incorrect format for tagged value\n";
                    std::cerr << "Field: " << str << "\n";
                    exit(EXIT_FAILURE);
                }

                if(tokens[1].size() != 1 && tokens[1][0] != getTypeCode(m_value))
                {
                    std::cerr << "Unexpected type in tagged value. Expected: " << getTypeCode(m_value);
                    std::cerr << " Found: " << tokens[1] << "\n";
                    exit(EXIT_FAILURE);
                }

                std::stringstream ss(tokens[2]);
                ss >> m_value;
                m_isInitialized = true;
            }
        
        private:
            T m_value;
            bool m_isInitialized;
    };

    // These are the valid tags that can be used
    typedef TagValue<char> CharTag;
    typedef TagValue<int> IntTag;
    typedef TagValue<float> FloatTag;
    typedef TagValue<std::string> StringTag;
};

#endif
