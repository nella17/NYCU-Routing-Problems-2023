#include "utils.hpp"

std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

int sign(int x) {
    return x == 0 ? 0 : x > 0 ? 1 : -1;
}

double sec_since(std::chrono::steady_clock::time_point start) {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count();
}
