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

#include "utils.hpp"

using namespace ISPDParser;

class GlobalRouter {
public:
    static const size_t C_SIZE = 11;

    using FP = void (GlobalRouter::*)(TwoPin*);

    struct Edge {
        int cap, demand, he, of;
        bool overflow;
        std::unordered_map<int, size_t> net;
        std::unordered_set<TwoPin*> twopins;
        Edge(int);
        void push(TwoPin*, int, int);
        void pop (TwoPin*, int, int);
    };

    struct Box {
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
    };

    bool stop;
    ISPDParser::ispdData* const ispdData;
    const std::array<ld, C_SIZE> C;

    GlobalRouter(ISPDParser::ispdData*, std::array<ld, C_SIZE>);

    void route();
    LayerAssignment::Graph* layer_assignment();

private:
    int k;
    size_t width, height;
    int min_width, min_spacing;
    std::vector<Edge> vedges, hedges;
    std::unordered_map<TwoPin*, Box> boxs;
    std::vector<TwoPin*> twopins;

    RPoint make(Point, Point);
    ld cost(TwoPin*);
    ld cost(Point, Point);
    ld cost(RPoint);
    ld cost(int, int, bool);
    ld cost(const Edge&) const;
    ld score(const TwoPin&) const;
    int delta(const TwoPin&) const;

    Edge& getEdge(RPoint);
    Edge& getEdge(int, int, bool);

    void ripup(TwoPin*);
    void place(TwoPin*);

    void Lshape(TwoPin*);
    void Lshape_impl(Path&, Point, Point);

    void Zshape(TwoPin*);

    void monotonic(TwoPin*);

    void VMR_impl(Point, Point, BoxCost&);
    void HMR_impl(Point, Point, BoxCost&);
    void HUM_impl(Path&, Point, Point, Box);
    void HUM(TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void init();
    int check_overflow();

    void ripup_place(FP);
    void routing(const char*, FP, int = 1);

    void print_edges();
};
