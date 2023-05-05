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

using namespace ISPDParser;

class GlobalRouter {
public:
    static const size_t C_SIZE = 10;

    ISPDParser::ispdData* const ispdData;
    const std::array<ld, C_SIZE> C;
    GlobalRouter(ISPDParser::ispdData*, std::array<ld, C_SIZE>);

    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    struct Edge {
        int cap, demand, he, of;
        std::map<int, size_t> net;
        std::set<TwoPin*> twopins;
        Edge(int);
        void push(TwoPin*, int, int);
        void pop (TwoPin*, int, int);
    };

    struct Box {
        int L, R, U, D;
        std::array<Point,4> points() const;
    };

    int k;
    ld cost(const Edge&);
    ld score(const TwoPin&);
    int delta(const TwoPin&);

    size_t width, height;
    int min_width, min_spacing;
    std::vector<Edge> vedges, hedges;
    Edge& getEdge(int, int, bool);

    std::vector<TwoPin*> twopins;
    void ripup(TwoPin*);
    void place(TwoPin*);

    Path Lshape(Point, Point);
    Path HUM(Point, Point);

    void construct_2D_grid_graph();
    void net_decomposition();
    void init_congestion();
    void pattern_routing();
    int HUM_routing();
};
