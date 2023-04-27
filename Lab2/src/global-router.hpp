#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

class GlobalRouter {
public:
    ISPDParser::ispdData* const ispdData;
    GlobalRouter(ISPDParser::ispdData*);
    void route(int);
    LayerAssignment::Graph* layer_assignment();

private:
    void construct_2D_grid_graph();
    void net_decomposition();
};
