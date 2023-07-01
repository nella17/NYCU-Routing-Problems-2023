#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include "src/ispdData.h"
#include "src/grid-graph.hpp"
#define _ <<' ' <<
#define ALL(v) v.begin(),v.end()

struct Edge {
    int cap, demand, via;
    Edge(int _c = 0): cap(_c), demand(0), via(0) {}
};

template<typename T>
class Grid {
    size_t w, h;
    std::vector<T> nodes;

public:
    auto width()  { return w; }
    auto height() { return h; }

    inline size_t rp2idx(int x, int y) const {
        return (size_t)x * h + (size_t)y;
    }

    const T& at(int x, int y) const {
        return nodes.at(rp2idx(x, y));
    }

    T& at(int x, int y) {
        return nodes.at(rp2idx(x, y));
    }

    void init(size_t width, size_t height, T I = T()) {
        w = width; h = height;
        nodes.clear();
        nodes.assign(w * h, I);
    }

    auto begin() { return nodes.begin(); }
    auto begin() const { return nodes.begin(); }
    auto end() { return nodes.end(); }
    auto end() const { return nodes.end(); }
};

signed main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "%s: <input> <output> <cmap>", argv[0]);
        return -1;
    }

    auto in   = argv[1];
    auto out  = argv[2];
    auto cmap = argv[3];

    // parse input
    std::ifstream fin(in);
    if (!fin.is_open()) {
        fprintf(stderr, "Failed to open input file\n");
        exit(EXIT_FAILURE);
    }
    auto ispdData = ISPDParser::parse(fin);
    fin.close();
    // std::cout << *ispdData << std::endl;

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ALL(ispdData->nets), [&](auto net) {
        for (auto &_pin : net->pins) {
            int x = (std::get<0>(_pin) - ispdData->lowerLeftX) / ispdData->tileWidth;
            int y = (std::get<1>(_pin) - ispdData->lowerLeftY) / ispdData->tileHeight;
            int z = std::get<2>(_pin) - 1;
            if (std::any_of(ALL(net->pin3D), [x, y, z](const auto &pin) {
                return pin.x == x && pin.y == y && pin.z == z;
            })) continue;
            net->pin3D.emplace_back(x, y, z);
            if (std::any_of(ALL(net->pin2D), [x, y](const auto &pin) { 
                return pin.x == x && pin.y == y;
            })) continue;
            net->pin2D.emplace_back(x, y);
        }
        return net->pin3D.size() > 1000 || net->pin2D.size() <= 1;
    }), ispdData->nets.end());
    ispdData->numNet = (int)ispdData->nets.size();

    auto width  = (size_t)ispdData->numXGrid;
    auto height = (size_t)ispdData->numYGrid;
    auto verticalCapacity   = std::accumulate(ALL(ispdData->verticalCapacity), 0);
    auto horizontalCapacity = std::accumulate(ALL(ispdData->horizontalCapacity), 0);

    GridGraph<Edge> gridedge;
    gridedge.init(width, height, Edge(verticalCapacity), Edge(horizontalCapacity));
    Grid<int> grid;
    grid.init(width, height, 0);
    for (auto capacityAdj: ispdData->capacityAdjs) {
        auto [x1,y1,z1] = capacityAdj->grid1;
        auto [x2,y2,z2] = capacityAdj->grid2;
        if (z1 != z2) continue;
        auto z = (size_t)z1 - 1;
        auto lx = std::min(x1, x2), rx = std::max(x1, x2);
        auto ly = std::min(y1, y2), ry = std::max(y1, y2);
        auto dx = rx - lx, dy = ry - ly;
        if (dx + dy != 1) continue;
        auto hori = dx;
        auto& edge = gridedge.at(lx, ly, hori);
        auto layerCap = (hori ? ispdData->horizontalCapacity : ispdData->verticalCapacity)[z];
        edge.cap -= layerCap - capacityAdj->reducedCapacityLevel;
        // std::cerr _ dx _ dy _ "/" _ lx _ ly _ "/" _ cong.cap _ layerCap _ std::endl;
    }

    // parse output
    std::ifstream fout(out);
    if (!fout.is_open()) {
        fprintf(stderr, "Failed to open output file\n");
        exit(EXIT_FAILURE);
    }
    while (!fout.eof()) {
        std::string netname, line;
        int id;
        fout >> netname >> id; fout.ignore();
        while (getline(fout, line) and line != "!") {
            // std::cerr _ line _ std::endl;
            std::stringstream ss(line);
            char c;
            int sx, sy, sz, ex, ey, ez;
            ss
                >> c >> sx >> c >> sy >> c >> sz >> c
                >> c
                >> c >> ex >> c >> ey >> c >> ez >> c
                ;

            sx = (sx - ispdData->lowerLeftX) / ispdData->tileWidth;
            sy = (sy - ispdData->lowerLeftY) / ispdData->tileHeight;
            sz = sz - 1;

            ex = (ex - ispdData->lowerLeftX) / ispdData->tileWidth;
            ey = (ey - ispdData->lowerLeftY) / ispdData->tileHeight;
            ez = ez - 1;

            if (sx > ex) std::swap(sx, ex);
            if (sy > ey) std::swap(sy, ey);
            if (sz > ez) std::swap(sz, ez);

            if (sz == ez) {
                // edge
                if (sx == ex) {
                    auto x = sx;
                    for (int y = sy; y < ey; y++)
                        gridedge.at(x, y, 0).demand++;
                } else {
                    auto y = sy;
                    for (int x = sx; x < ex; x++)
                        gridedge.at(x, y, 1).demand++;
                }
            } else {
                // via
                grid.at(sx, sy)++;
            }
        }
    }
    fout.close();

    auto L = 0, R = (int)width-1, B = 0, U = (int)height-1;
    /*
    std::cerr << "horizontal edges\n";
    for (int j = U; j >= B; j--) {
        for (int i = L; i+1 <= R; i++) {
            auto& e = gridedge.at(i, j, 1);
            std::cerr _ '(' << i << ',' << j << ')' << e.demand <<'/'<< e.cap;
        }
        std::cerr _ std::endl;
    }
    std::cerr << "vertical edges\n";
    for (int j = U-1; j >= B; j--) {
        for (int i = L; i <= R; i++) {
            auto& e = gridedge.at(i, j, 0);
            std::cerr _ '(' << i << ',' << j << ')' << e.demand <<'/'<< e.cap;
        }
        std::cerr _ std::endl;
    }
    //*/

    Grid<Edge> image;
    auto iw = width * 2 - 1, ih = height * 2 - 1;
    image.init(ih, iw);

    for (size_t i = 0; i < ih; i++) for (size_t j = 0; j < iw; j++) {
        auto x = int(j / 2), y = int(i / 2);
        auto& data = image.at((int)i, (int)j);
        switch (((i % 2) * 2) + (j % 2)) {
            case 0: // node
                data.via = grid.at(x, y);
                break;
            case 1: // h edge
                data = gridedge.at(x, y, 1);
                data.via = -1;
                break;
            case 2: // v edge
                data = gridedge.at(x, y, 0);
                data.via = -2;
                break;
            default: // other
                data.via = -3;
        }
    }

    std::ofstream fcmap(cmap);
    if (!fcmap.is_open()) {
        fprintf(stderr, "Failed to open input file\n");
        exit(EXIT_FAILURE);
    }
    fcmap << ispdData->numLayer << '\n';
    for (size_t i = 0; i < ih; i++) for (size_t j = 0; j < iw; j++) {
        auto& e = image.at((int)(ih-1-i), (int)j);
        fcmap << e.demand << '/' << e.cap << '/' << e.via << " \n"[j+1==iw];
    }
    fcmap.close();

    return 0;
}
