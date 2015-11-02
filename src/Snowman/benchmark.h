#ifndef SNOWMAN_BENCHMARK_H__
#define SNOWMAN_BENCHMARK_H__

#include <string>
#include <vector>

void parseBenchmarkOptions(int argc, char** argv);
void runBenchmark(int argc, char** argv);
void assemblyTest();
void sampleReads(const std::string& seq, std::vector<std::string>& reads, int cov, double error_rate);
void realignRandomSegments();
std::string genBreaks();
std::vector<double> parseErrorRates(const std::string& s);
std::string errorRateString(const std::vector<double>& v, const std::string& name);
void splitBam();

#endif
