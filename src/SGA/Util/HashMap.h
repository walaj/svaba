//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HashMap compability defines. Similar to Shaun Jackman's
// compatability code in abyss (see Common/HashMap)
//
//
#ifndef HASHMAP_H
#define HASHMAP_H

#include "config.h"

#if HAVE_UNORDERED_MAP

# include <unordered_map>
# include <unordered_set>
#define HashMap std::unordered_map
#define HashSet std::unordered_set
typedef std::hash<std::string> StringHasher;
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#define HashMap std::tr1::unordered_map
#define HashSet std::tr1::unordered_set
typedef std::tr1::hash<std::string> StringHasher;
#elif HAVE_EXT_HASH_MAP
#define USING_EXT_HASH_MAP 1
# undef __DEPRECATED
#include <ext/hash_map>
#include <ext/hash_set>
#define HashMap __gnu_cxx::hash_map
#define HashSet __gnu_cxx::hash_set
#else
# error No hash map implementation found
#endif

# if USING_EXT_HASH_MAP
# include <cstddef>
# include <string>
struct StringHasher
{
    size_t operator()( const std::string& x) const
    {
        return __gnu_cxx::hash<const char*>()( x.c_str());
    }
};
#endif 

// Ensure the sparse hash is available
#if HAVE_GOOGLE_SPARSE_HASH_MAP
#include <google/sparse_hash_map>
//#define SparseHashMap google::sparse_hash_map // JEREMIAH
#else
#define SparseHashMap std::tr1::unordered_map
//#error The google sparse hash is required
#endif 

#endif
