#pragma once

// #undef DEBUG
#define DEBUG 1
#if DEBUG
#include <cassert>
#else
#define assert(x) (x)
#endif

#define _ <<' '<<
#define ALL(v) v.begin(), v.end()
