#include "global-router.hpp"

#include <queue>

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData): ispdData(_ispdData) {}

void GlobalRouter::route(int timeLimitSec) {
    construct_2D_grid_graph();
    net_decomposition();
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
