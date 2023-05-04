#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include <vector>
#include <map>

class GlobalRouter {
public:
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

    struct CongSnapshot {
        size_t width, height, min_width, min_spacing;
        std::vector<EdgeInfo> vcong, hcong;
    };
    struct Congestion {
        size_t width, height, min_width, min_spacing;
        std::vector<Edge> vcong, hcong;
        CongSnapshot dump() const;
        void init(ISPDParser::ispdData*);
        Edge& getEdge(int, int, bool);
    };

    ISPDParser::ispdData* const ispdData;
    GlobalRouter(ISPDParser::ispdData*);

    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    Congestion congestion;

    std::vector<ISPDParser::TwoPin*> twopins;
    void ripup(ISPDParser::TwoPin*);
    void place(ISPDParser::TwoPin*);

    void Lshape(ISPDParser::TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void init_route();
};
