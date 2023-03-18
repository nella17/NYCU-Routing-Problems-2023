#pragma once

#if DEBUG
#include <cassert>
#else
#define assert(x) ((x) ? 0 : throw 1)
#endif

#define _ <<' '<<
#define ALL(v) v.begin(), v.end()
