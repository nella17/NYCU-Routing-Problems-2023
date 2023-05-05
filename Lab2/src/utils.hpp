#pragma once

#include <iostream>
#include <random>
#include <chrono>

#define _ <<' ' <<
#define ALL(v) v.begin(),v.end()

using ld = long double;

extern std::mt19937 rng;

int sign(int);

double sec_since(std::chrono::steady_clock::time_point);

template<typename T>
T randint(T l, T r) { return std::uniform_int_distribution<T>(l,r)(rng); }

template<typename T>
T randint(T n) { return randint(0,n-1); }

template<typename T>
inline T average(const std::vector<T>& v) {
    T i{};
    for (const auto& x: v)
        i += x;
    return i / (T)v.size();
}
