#pragma once

#include "ispdData.h"
#include "LayerAssignment.h"

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
