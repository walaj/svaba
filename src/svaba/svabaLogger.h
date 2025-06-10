#pragma once
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

class SvabaOptions;

/// A simple thread safe logger that writes to a file and (optionally) stderr.
class SvabaLogger {
public:
  SvabaLogger() = default;
  ~SvabaLogger();

  /// Must be called once at program startup, before any log() calls.
  /// Opens (and appends to) the given filename.
  void init(const std::string& filename);

  /// Log a message.  
  /// @param toErr if true, also write to std::cerr  
  /// @param toLog if true, write to the log file  
  /// @param args  one or more things you can stream into an ostringstream
  template <typename... Args>
  void log(bool toErr, bool toLog, Args&&... args) {

    // it no printing at all, just leave
    if (!toErr && !toLog)
      return;

    std::lock_guard<std::mutex> lock(mtx_);
    std::ostringstream oss;
    // fold expression to stream all args
    (oss << ... << std::forward<Args>(args));
    auto msg = oss.str();
    if (toLog && logFile_.is_open()) logFile_ << msg << "\n";
    if (toErr) std::cerr << msg << "\n";
  }

  void welcome(SvabaOptions& opts);
  
private:
  std::ofstream   logFile_;
  std::mutex      mtx_;
};
