#include "utils.hpp"

std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());
