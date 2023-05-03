#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include <vector>

class GlobalRouter {
public:
    ISPDParser::ispdData* const ispdData;
    GlobalRouter(ISPDParser::ispdData*);
    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    size_t width, height;

    std::vector<ISPDParser::TwoPin*> twopins;
    void init(ISPDParser::TwoPin*);
    void add(ISPDParser::TwoPin*);

    void construct_2D_grid_graph();
    void net_decomposition();
};
