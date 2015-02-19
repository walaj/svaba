#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include "GenomicRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <memory>

typedef std::shared_ptr<BamTools::BamAlignment> BamAlignmentUP;

namespace SnowUtils {

template <typename T> 
std::string AddCommas(T data) { 
  std::stringstream ss; ss << data; std::string s = ss.str();
  if (s.length() > 3)
     for (int i = s.length()-3; i > 0; i -= 3)
       s.insert(i,",");
   return s;
}

using std::ifstream;
inline bool existTest (const std::string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}

// clean an output directory
/*inline std::string getDirPath(std::string dir) {
 
  // get the basepath from directory
  // check if it is a directory
  struct stat sb;
  std::string odir = dir;
  if (stat(dir.c_str(), &sb) == 0 && !S_ISDIR(sb.st_mode)) // it is successful, but not directory
     odir  = dir.substr(0, dir.find_last_of('/'));

  return odir;
 
  }*/

inline void displayRuntime(const timespec start) {

  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double elapsed = (finish.tv_sec - start.tv_sec);
  int t = clock()/CLOCKS_PER_SEC;
  int min = (int)floor(elapsed / 60.0);
  int sec = (int)(elapsed-min*60);
  char buffer[100];
  sprintf (buffer, "CPU: %4dm%02ds Wall: %4dm%02ds", 
            (int)floor( ((double)t) /60.0), t % 60, min, sec);
  printf ("%s",buffer);
}

inline void rcomplement(std::string &a) {

  std::reverse(&a[0], &a[a.size()]);
  std::string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';
}

  // calculate the percentage
 template <typename T> inline int percentCalc(T numer, T denom) {
   if (denom <= 0)
     return 0;
   int perc  = static_cast<int>(floor((float)numer / (float)denom * 100.0));
   return perc;
 }

 // remove the last character from a string
 inline std::string cutLastChar(std::string in) {
   if (in.length() == 0)
     return in;
   else 
     return in.substr(0, in.length() - 1);
 }

 // remove substrings from a string
 inline std::string scrubString(std::string toscrub, std::string toremove) {
   std::string::size_type i = toscrub.find(toremove);
   while (i != std::string::npos) {
     toscrub.erase(i, toremove.length());
     i = toscrub.find(toremove);
   }
   return toscrub;
 }

 // get a file name + extension
 // https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch10s16.html
 inline std::string getFileName(const std::string& s) {

   char sep = '/';

   #ifdef _WIN32
     sep = '\\';
   #endif

   size_t i = s.rfind(sep, s.length());
   if (i != std::string::npos && ( (i+1) < s.length()) ) {
     return(s.substr(i+1, s.length()));
   }

   return(s);
 }

 // get a string or int tag that might be separted by "x"
 std::vector<std::string> GetStringTag(const BamAlignmentUP& a, const std::string tag);
 std::vector<int> GetIntTag(const BamAlignmentUP& a, const std::string tag);

 // add a tag, and if its already there separate by "x"
 void SmartAddTag(BamAlignmentUP &a, const std::string tag, const std::string val);

} // end namespace

#endif
