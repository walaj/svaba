#include "svabaLogger.h"
#include "svaba_params.h"

#include "svabaOptions.h"

SvabaLogger::~SvabaLogger() {
  if (logFile_.is_open()) logFile_.close();
}

void SvabaLogger::init(const std::string& filename) {
  logFile_.open(filename, std::ios::out | std::ios::app);
  if (!logFile_) {
    std::cerr << "ERROR: Unable to open log file: " << filename << std::endl;
  }
}

void SvabaLogger::welcome(SvabaOptions& opts) {
  // write to log
  this->log(true, true, 
	    "-----------------------------------------------------------\n",
	    "---  Running svaba SV and indel detection on ",
	    " threads ---", (opts.numThreads >= 10 ? "\n" : "-\n"),
	    "---  Version: ", SVABA_VERSION, " - ", SVABA_DATE, "                           ---\n",
	    "---    (inspect *.log for real-time progress updates)   ---\n",
	    "-----------------------------------------------------------");
  
  // print the options to console if verbose (handled insize SvabaOptions)
  //opts.printToLogger(*this);
}
