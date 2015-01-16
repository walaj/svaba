#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>

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

} // end namespace

#endif
