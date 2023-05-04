#include "global-router.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>
#include <random>
#include <chrono>

#include <iostream>
#define _ <<' ' <<
#define ALL(v) v.begin(),v.end()

std::mt19937 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());
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

GlobalRouter::Edge::Edge(int _cap): cap(_cap), net{}, twopins{} {}

GlobalRouter::EdgeInfo GlobalRouter::Edge::dump() const  {
    return EdgeInfo{
        .cap = this->cap,
        .daemon = (int)this->net.size(),
    };
}

GlobalRouter::Edge& GlobalRouter::getEdge(int x, int y, bool hori) {
    if (hori)
        return hcong[ (size_t)x * height + (size_t)y ];
    else
        return vcong[ (size_t)x + (size_t)y * width ];
}

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData): ispdData(_ispdData) {}

void GlobalRouter::ripup(ISPDParser::TwoPin* twopin) {
    assert(!twopin->ripup);
    // TODO
}

void GlobalRouter::place(ISPDParser::TwoPin* twopin) {
    assert(twopin->ripup);
    // TODO
}

void GlobalRouter::Lshape(ISPDParser::TwoPin* twopin) {
    assert(twopin->ripup);
    twopin->path.clear();
    auto f = twopin->from, t = twopin->to;
    if (f.y > t.y) std::swap(f, t);
    if (f.x > t.x) std::swap(f, t);
    auto Lx = [&](int y, int lx, int rx) {
        if (lx > rx) std::swap(lx, rx);
        for (auto x = lx; x < rx; x++)
            twopin->path.emplace_back(x, y, 1);
    };
    auto Ly = [&](int x, int ly, int ry) {
        if (ly > ry) std::swap(ly, ry);
        for (auto y = ly; y < ry; y++)
            twopin->path.emplace_back(x, y, 0);
    };
    auto L = [&](ISPDParser::Point p1, ISPDParser::Point p2) {
        // std::cerr _ "L" _ p1 _ p2 _ std::endl;
        if (p1.x == p2.x) return Ly(p1.x, p1.y, p2.y);
        if (p1.y == p2.y) return Lx(p1.y, p1.x, p2.x);
    };
    int d = randint(2);
    ISPDParser::Point m(
        d ? f.x : t.x,
        d ? t.y : f.y
    );
    L(f, m); L(m, t);
}

void GlobalRouter::route(int timeLimitSec) {
    construct_2D_grid_graph();
    net_decomposition();
    for (auto net: ispdData->nets)
        for (auto& twopin: net->twopin)
            twopins.emplace_back(&twopin);
    init_route();
    // exit(-1);
}

void GlobalRouter::construct_2D_grid_graph() {
    // construct 2D grid graph

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ALL(ispdData->nets), [&](ISPDParser::Net *net) {

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
}

void GlobalRouter::net_decomposition() {
    for (auto net: ispdData->nets) {
        auto sz = net->pin2D.size();
        net->twopin.clear();
        net->twopin.reserve(sz - 1);
        std::vector<bool> vis(sz, false);
        std::priority_queue<std::tuple<int, size_t, size_t>> pq{};
        auto add = [&](size_t i) {
            vis[i] = true;
            auto [xi, yi, zi] = net->pin2D[i];
            for (size_t j = 0; j < sz; j++) if (!vis[j]) {
                auto [xj, yj, zj] = net->pin2D[j];
                auto d = std::abs(xi - xj) + std::abs(yi - yj);
                pq.emplace(-d, i, j);
            }
        };
        add(0);
        while (!pq.empty()) {
            auto [d, i, j] = pq.top(); pq.pop();
            if (vis[j]) continue;
            net->twopin.emplace_back(
                net->pin2D[i],
                net->pin2D[j],
                net
            );
            add(j);
        }
        assert(net->twopin.size() == sz - 1);
        assert(net->twopin.capacity() == sz - 1);
    }
}

void GlobalRouter::init_congestion() {
    width  = (size_t)ispdData->numXGrid;
    height = (size_t)ispdData->numYGrid;
    min_width = (size_t)average(ispdData->minimumWidth);
    min_spacing = (size_t)average(ispdData->minimumSpacing);
    auto verticalCapacity = std::accumulate(ALL(ispdData->verticalCapacity), 0);
    auto horizontalCapacity = std::accumulate(ALL(ispdData->horizontalCapacity), 0);
    vcong.assign(width * (height - 1), Edge(verticalCapacity));
    hcong.assign((width - 1) * height, Edge(horizontalCapacity));
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
        auto& cong = getEdge(lx, ly, hori);
        auto layerCap = (hori ? ispdData->horizontalCapacity : ispdData->verticalCapacity)[z];
        cong.cap -= layerCap - capacityAdj->reducedCapacityLevel;
        // std::cerr _ dx _ dy _ "/" _ lx _ ly _ "/" _ cong.cap _ layerCap _ std::endl;
    }
    /*
    std::cerr << "horizontalCapacity\n";
    for (int j = height-1; j >= 0; j--)
        for (int i = 0; i+1 < width; i++)
            std::cerr << getCong(i, j, 1).cap << " \n"[i+2==width];
    std::cerr << "verticalCapacity\n";
    for (int j = height-2; j >= 0; j--)
        for (int i = 0; i < width; i++)
            std::cerr << getCong(i, j, 0).cap << " \n"[i+1==width];
    //*/
}

void GlobalRouter::init_route() {
    // auto savedCong
    for (auto twopin: twopins) {
        twopin->ripup = true;
        twopin->reroute = false;
        Lshape(twopin);
        // std::cerr _ *twopin _ std::endl;
        // for (auto x
        place(twopin);
    }
}

LayerAssignment::Graph* GlobalRouter::layer_assignment() {
    // Assign routing layers to the two-pin net
    auto graph = new LayerAssignment::Graph;
    graph->initialLA(*ispdData, 1);
    graph->convertGRtoLA(*ispdData, true);
    graph->COLA(true);
    return graph;
}
