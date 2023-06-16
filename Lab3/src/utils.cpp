#include "utils.hpp"

#include <algorithm>

std::vector<std::vector<bool>> comb(size_t n, size_t k) {
    std::vector<std::vector<bool>> r{};
    if (k == 0 or k > n) return r;
    std::vector<bool> bitmask(n, false);
    for (size_t i = 0; i < k; i++)
        bitmask[i] = true;
    size_t size = 1;
    for (size_t i = 0; i < k; i++) {
        size *= n - i;
        size /= i+1;
    }
    r.reserve(size);
    do {
        r.emplace_back(bitmask);
    } while (std::prev_permutation(ALL(bitmask)));
    return r;
}

std::vector<std::vector<size_t>> comb_id(size_t n, size_t k) {
    auto c = comb(n, k);
    std::vector<std::vector<size_t>> r{};
    r.reserve(c.size());
    for (auto&& mask: c) {
        std::vector<size_t> v{};
        v.reserve(k);
        for (size_t i = 0; i < n; i++)
            if (mask[i])
                v.emplace_back(i);
        r.emplace_back(v);
    }
    return r;
}
