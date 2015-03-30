#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include "GenomicRegion.h"
#include <unistd.h>

#include "reads.h"

typedef std::vector<BamTools::CigarOp> CigarOpVec;

namespace SnowUtils {

inline bool read_access_test (const std::string& name) {
  return (access (name.c_str(), R_OK) == 0); 
}

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
/* inline std::string getDirPath(std::string dir) {
 
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
 /* inline std::string getFileName(const std::string& s) {

   char sep = '/';

   #ifdef _WIN32
     sep = '\\';
   #endif

   size_t i = s.rfind(sep, s.length());
   if (i != std::string::npos && ( (i+1) < s.length()) ) {
     return(s.substr(i+1, s.length()));
   }

   return(s);
   }*/

 /*! @function parse a tag storing multiple strings separated by 'x' character
  * @param Read containing tag to be parsed
  * @param tag to parse
  * @return all of the strings from the tag
  */
 std::vector<std::string> GetStringTag(const Read& a, const std::string tag);

 /*! @function parse a tag storing multiple integers separated by 'x' character
  * @param Read containing tag to be parsed
  * @param tag to parse
  * @return all of the integers from the tag
  */
 std::vector<int> GetIntTag(const Read& a, const std::string tag);

 // add a tag, and if its already there separate by "x"
 void SmartAddTag(Read &a, const std::string tag, const std::string val);

 /*! @function Convert a CigarOpVec to a string for printing
  * @param Cigar vec to be read
  * @return CIGAR string
  */
 std::string cigarToString(const CigarOpVec &cig);

 /*! @function Flip the cigar so that is in opposite orientation
  * @param cigar to be flipped in place
  */
 void flipCigar(CigarOpVec &cig);

 /*! @function Parse a cigar string into a vector<CigarOp>
  * @param CIGAR string to be parsed
  * @return parsed CIGAR in vector format from BamTools package
  */
 CigarOpVec stringToCigar(const std::string& val);

 /*! @function Parse tags from a SAM alignment and add to a BamAlignment
  * @param tag to be parsed (e.g. XA:Z:...)
  * @param alignment object to be modified
  */
 void parseTags(const std::string& val, BamTools::BamAlignment &a);

}

#endif
