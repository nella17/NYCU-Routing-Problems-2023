#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include <vector>
#include <unordered_set>

class GlobalRouter {
public:
    struct Congestion {
        int cap;
        std::multiset<size_t> net{};
        Congestion(int);
    };

    ISPDParser::ispdData* const ispdData;
    GlobalRouter(ISPDParser::ispdData*);

    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    size_t width, height;

    std::vector<Congestion> vcong, hcong;
    Congestion& getCong(int, int, bool);

    std::vector<ISPDParser::TwoPin*> twopins;
    void init(ISPDParser::TwoPin*);
    void add(ISPDParser::TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
    void init_cap();
};
