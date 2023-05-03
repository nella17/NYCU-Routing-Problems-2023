#include "global-router.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>

#include <iostream>
#define _ <<' ' <<

GlobalRouter::Congestion::Congestion(int _cap): cap(_cap), net{} {}

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData): ispdData(_ispdData) {}

GlobalRouter::Congestion& GlobalRouter::getCong(int x, int y, bool hori) {
    if (hori)
        return hcong[ (size_t)x * height + (size_t)y ];
    else
        return vcong[ (size_t)x + (size_t)y * width ];
}

void GlobalRouter::init(ISPDParser::TwoPin* twopin) {
    // TODO
}

void GlobalRouter::add(ISPDParser::TwoPin* twopin) {
    // TODO
}

void GlobalRouter::route(int timeLimitSec) {
    width  = (size_t)ispdData->numXGrid;
    height = (size_t)ispdData->numYGrid;
    construct_2D_grid_graph();
    net_decomposition();
    for (auto net: ispdData->nets)
        for (auto& twopin: net->twopin)
            twopins.emplace_back(&twopin);
    init_cap();
}

void GlobalRouter::construct_2D_grid_graph() {
    // construct 2D grid graph

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ispdData->nets.begin(), ispdData->nets.end(), [&](ISPDParser::Net *net) {

        for (auto &_pin : net->pins) {

            int x = (std::get<0>(_pin) - ispdData->lowerLeftX) / ispdData->tileWidth;
            int y = (std::get<1>(_pin) - ispdData->lowerLeftY) / ispdData->tileHeight;
            int z = std::get<2>(_pin) - 1;

            if (std::any_of(net->pin3D.begin(), net->pin3D.end(), [x, y, z](const auto &pin) {
                return pin.x == x && pin.y == y && pin.z == z;
            })) continue;
            net->pin3D.emplace_back(x, y, z);

            if (std::any_of(net->pin2D.begin(), net->pin2D.end(), [x, y](const auto &pin) { 
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

void GlobalRouter::init_cap() {
    auto verticalCapacity = std::accumulate(
        ispdData->verticalCapacity.begin(),
        ispdData->verticalCapacity.end(),
        0
    );
    auto horizontalCapacity = std::accumulate(
        ispdData->horizontalCapacity.begin(),
        ispdData->horizontalCapacity.end(),
        0
    );
    vcong.assign(width * (height - 1), Congestion(verticalCapacity));
    hcong.assign((width - 1) * height, Congestion(horizontalCapacity));
    for (auto capacityAdj: ispdData->capacityAdjs) {
        auto [x1,y1,z1] = capacityAdj->grid1;
        auto [x2,y2,z2] = capacityAdj->grid2;
        assert(z1 == z2);
        auto z = (size_t)z1 - 1;
        auto lx = std::min(x1, x2), rx = std::max(x1, x2);
        auto ly = std::min(y1, y2), ry = std::max(y1, y2);
        auto dx = rx - lx, dy = ry - ly;
        assert(dx + dy == 1);
        auto hori = dx;
        auto& cong = getCong(lx, ly, hori);
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

LayerAssignment::Graph* GlobalRouter::layer_assignment() {
    // Assign routing layers to the two-pin net
    auto graph = new LayerAssignment::Graph;
    graph->initialLA(*ispdData, 1);
    graph->convertGRtoLA(*ispdData, true);
    graph->COLA(true);
    return graph;
}
