#pragma once

// Re-run LOD / PASS filtering against an existing bps.txt.gz. See refilter.cpp
// for option parsing and the supported input layouts (41-col legacy and
// 51-col SvABA2.0). Invoked from svaba.cpp as `svaba refilter ...`.
void runRefilterBreakpoints(int argc, char** argv);
