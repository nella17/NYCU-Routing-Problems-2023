#pragma once

#include <vector>
#include <cstddef>
#define _ <<' '<<
#define ALL(v) v.begin(),v.end()

std::vector<std::vector<bool>> comb(size_t n, size_t k);
std::vector<std::vector<size_t>> comb_id(size_t n, size_t k);
