#include "snowUtil.h"
#include <cmath>
#include <cstdio>
#include "GenomicRegion.h"

using namespace std;
using namespace BamTools;

void displayRuntime(const timespec start) {
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

// get the SamHeader from a bamfile                                                                                                                                                                                 
string getSamHeader(string bamfile, SamHeader &sam) {                                                                                                                                                    
                                                                                                                                                                                                                    
  BamReader read;                                                                                                                                                                                                   
  if (!read.Open(bamfile)) {                                                                                                                                                                                        
    cerr << "Failed to open contig BAM to get header on bam " << bamfile << endl;                                                                                                                                   
    cerr << "   Setting default of 'none'" << endl;                                                                                                                                                                 
    sam = SamHeader("none");                                                                                                                                                                                        
    read.Close();                                                                                                                                                                                                   
  }                                                                                                                                                                                                                 
                                                                                                                                                                                                                    
  sam = read.GetHeader();                                                                                                                                                                                           

  return read.GetHeaderText();
                                                                                                                                                                                                                    
}                                                                                                                                                                                                                   

// get the reference vector from a BAM header                                                                                                                                                                       
void getRefVector(string bamfile, RefVector &ref) {                                                                                                                                                    
                                                                                                                                                                                                                    
  BamReader read;                                                                                                                                                                                                   
  if (!read.Open(bamfile)) {                                                                                                                                                                                        
    cerr << "Failed to open contig BAM to get reference on bam " << bamfile << endl;                                                                                                                                
    cerr << "   Setting default of 1,2,..., Y" << endl;                                                                                                                                                             
    for (int i = 0; i < 25; i++) {                                                                                                                                                                                  
      RefData rf(CHR_NAME[i], CHR_LEN[i]);                                                                                                                                                                          
      ref.push_back(rf);                                                                                                                                                                                            
    }                                                                                                                                                                                                               
    read.Close();                                                                                                                                                                                                   
    return;                                                                                                                                                                                                         
  }                                                                                                                                                                                                                 
                                                                                                                                                                                                                    
  ref = read.GetReferenceData();                                                                                                                                                                                    
  read.Close();                                                                                                                                                                                                     
} 
