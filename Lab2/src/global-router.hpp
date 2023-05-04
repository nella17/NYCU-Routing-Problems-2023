#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include <array>
#include <vector>
#include <map>

#include "utils.hpp"

class GlobalRouter {
public:
    static const size_t C_SIZE = 8;

    struct EdgeInfo {
        int cap, daemon;
    };
    struct Edge {
        int cap;
        std::map<size_t, size_t> net;
        std::set<ISPDParser::TwoPin*> twopins;
        Edge(int);
        EdgeInfo dump() const;
    };

    ISPDParser::ispdData* const ispdData;
    const std::array<ld, C_SIZE> C;
    GlobalRouter(ISPDParser::ispdData*, std::array<ld, C_SIZE>);

    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    size_t width, height, min_width, min_spacing;
    std::vector<Edge> vcong, hcong;
    Edge& getEdge(int, int, bool);

    std::vector<ISPDParser::TwoPin*> twopins;
    void ripup(ISPDParser::TwoPin*);
    void place(ISPDParser::TwoPin*);

    void Lshape(ISPDParser::TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void init_congestion();
    void init_route();
};
