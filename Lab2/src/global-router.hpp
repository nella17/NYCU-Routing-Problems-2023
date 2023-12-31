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
        int used; double cost;
        Edge(int);
        inline bool overflow() const;
        inline bool push(TwoPin*);
        inline bool pop (TwoPin*);
    };

    struct Box {
        bool eL, eR, eB, eU;
        int L, R, B, U;
        Box(Point, Point);
        inline Point BL() const;
        inline Point UR() const;
        inline size_t width() const;
        inline size_t height() const;
    };
    struct BoxCost: Box {
        struct Data {
            double cost = INFINITY;
            std::optional<Point> from = std::nullopt;
        };
        std::vector<Data> cost;
        BoxCost(const Box&);
        inline Data& operator()(Point);
        inline Data& operator()(int,int);
        inline void trace(Path&, Point);
    };

    struct Net {
        int overflow, overflow_twopin, wlen, reroute;
        double score, cost;
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
    int min_width, min_spacing, min_net, mx_cap;
    GridGraph<Edge> grid;
    std::vector<Net*> nets;
    std::vector<TwoPin*> twopins;
    // TODO
    // std::unordered_map<int, Net*> id2net;

    int selcost;
    inline double cost(const TwoPin*) const;
    inline double cost(Point, Point) const;
    inline double cost(RPoint) const;
    inline double cost(int, int, bool) const;
    inline double cost(const Edge&) const;

    void del_cost(Net*);
    void del_cost(TwoPin*);
    void add_cost(Net*);
    void add_cost(TwoPin*);
    double calc_cost(const Edge&) const;

    static constexpr int COSTSZ  = 1024;
    static constexpr int COSTOFF = 256;
    double cost_pe[COSTSZ];
    inline double get_cost_pe(int) const;
    void build_cost();

    void sort_twopins();
    inline double score(const TwoPin*) const;
    inline double score(const Net*) const;
    inline int delta(const TwoPin*) const;

    inline const Edge& getEdge(RPoint) const;
    inline Edge& getEdge(RPoint);
    inline const Edge& getEdge(int, int, bool) const;
    inline Edge& getEdge(int, int, bool);

    void ripup(TwoPin*);
    void place(TwoPin*);

    void Lshape(TwoPin*);
    void Zshape(TwoPin*);

    inline void calcX(BoxCost&, int, int, int);
    inline void calcY(BoxCost&, int, int, int);

    void monotonic(TwoPin*);

    void VMR_impl(Point, Point, BoxCost&);
    void HMR_impl(Point, Point, BoxCost&);
    void HUM(TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void preroute();
    int check_overflow();

    void ripup_place(FP);
    void routing(const char*, FP, int = 1);

    void print_edges(int = -1, int = -1, int = -1, int = -1);
};
