#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include <array>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <optional>

#include "grid-graph.hpp"
#include "utils.hpp"

using namespace ISPDParser;

RPoint make(Point, Point);

class GlobalRouter {
public:
    using FP = void (GlobalRouter::*)(TwoPin*);

    struct Edge {
        int cap, demand, he, of;
        std::unordered_map<int, size_t> net;
        std::unordered_set<TwoPin*> twopins;
        Edge(int);
        bool overflow(int netId = -1) const;
        bool push(TwoPin*, int, int);
        bool pop (TwoPin*, int, int);
    };

    struct Box {
        bool eL, eR, eB, eU;
        int L, R, B, U;
        Box(Point, Point);
        Point BL() const;
        Point UR() const;
        size_t width() const;
        size_t height() const;
    };
    struct BoxCost: Box {
        struct Data {
            ld cost = INFINITY;
            std::optional<Point> from = std::nullopt;
        };
        std::vector<Data> cost;
        BoxCost(const Box&);
        Data& operator()(Point);
        Data& operator()(int,int);
        void trace(Path&, Point);
    };

    struct Net {
        int overflow, overflow_twopin, wlen, reroute;
        long double score, cost;
        ISPDParser::Net* net;
        std::vector<TwoPin*> twopins;
        Net(ISPDParser::Net*);
    };

    bool stop, print;
    ISPDParser::ispdData* const ispdData;

    GlobalRouter(ISPDParser::ispdData*);
    ~GlobalRouter();

    void route(bool = false);
    LayerAssignment::Graph* layer_assignment(bool = true);

private:
    int k;
    size_t width, height;
    int min_width, min_spacing, mx_cap;
    GridGraph<Edge> grid;
    std::unordered_map<TwoPin*, Box> boxs;
    std::vector<Net*> nets;
    std::vector<TwoPin*> twopins;
    std::unordered_map<int, Net*> id2net;

    ld cost(const TwoPin*) const;
    ld cost(int, Point, Point);
    ld cost(int, RPoint) const;
    ld cost(int, int, int, bool);
    ld cost(int, const Edge&) const;
    ld score(const TwoPin*) const;
    ld score(const Net*) const;
    int delta(const TwoPin*) const;

    const Edge& getEdge(RPoint) const;
    Edge& getEdge(RPoint);
    const Edge& getEdge(int, int, bool) const;
    Edge& getEdge(int, int, bool);

    void ripup(TwoPin*);
    void place(TwoPin*);

    void Lshape(TwoPin*);
    void Zshape(TwoPin*);

    void calcX(int, BoxCost&, int, int, int);
    void calcY(int, BoxCost&, int, int, int);

    void monotonic(TwoPin*);

    void VMR_impl(int, Point, Point, BoxCost&);
    void HMR_impl(int, Point, Point, BoxCost&);
    void HUM(TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void preroute();
    int check_overflow();

    void ripup_place(FP);
    void routing(const char*, FP, int = 1);

    void print_edges(int = -1, int = -1, int = -1, int = -1, int = -1);
};
