#pragma once

#include <vector>

template<typename T>
class GridGraph {
    size_t w, h;
    size_t vsz, hsz;
    std::vector<T> edges;

public:
    auto width()  { return w; }
    auto height() { return h; }

    inline size_t rp2idx(int x, int y, bool hori) const {
        if (hori)
            return (size_t)x * h + (size_t)y + vsz;
        else
            return (size_t)x + (size_t)y * w;
    }

    const T& at(int x, int y, bool hori) const {
        return edges.at(rp2idx(x, y, hori));
    }

    T& at(int x, int y, bool hori) {
        return edges.at(rp2idx(x, y, hori));
    }

    void init(size_t width, size_t height, T vI, T hI) {
        w = width; h = height;
        vsz = w * (h-1);
        hsz = (w-1) * h;
        edges.clear();
        edges.reserve(vsz + hsz);
        edges.insert(edges.end(), vsz, vI);
        edges.insert(edges.end(), hsz, hI);
    }

    auto begin() { return edges.begin(); }
    auto begin() const { return edges.begin(); }
    auto end() { return edges.end(); }
    auto end() const { return edges.end(); }
};
