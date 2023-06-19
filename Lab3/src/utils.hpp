#pragma once

#include <vector>
#include <queue>
#include <cstddef>
#define _ <<' '<<
#define ALL(v) v.begin(),v.end()

template<typename T>
using priority_queue_greater = std::priority_queue<T, std::vector<T>, std::greater<T>>;

std::vector<std::vector<bool>> comb(size_t n, size_t k);
std::vector<std::vector<size_t>> comb_id(size_t n, size_t k);
