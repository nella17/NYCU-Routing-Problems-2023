#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include "global-router.hpp"

#include <iostream>
#include <fstream>
#include <string>

auto parse_input(const char* file) {
    // parse input
    std::ifstream fp(file);
    if (!fp.is_open()) {
        fprintf(stderr, "Failed to open input file\n");
        exit(EXIT_FAILURE);
    }

    auto ispdData = ISPDParser::parse(fp);

    fp.close();
    // std::cout << *ispdData << std::endl;

    return ispdData;
}


int main(int argc, char* const argv []) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <inputFile> <outputFile> <timeLimitSec>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    auto inputFile = argv[1];
    auto outputFile = argv[2];
    auto timeLimitSec = argc < 4 ? 60 * 30 : atoi(argv[3]);

    auto ispdData = parse_input(inputFile);

    GlobalRouter gr(
        ispdData,
        { 7, 4, 150, 0.3, 30, 200, 30, 1 }
    );
    gr.route(timeLimitSec);

    // Describe the usage of the given layer assignment algorithm
    // Only works for the given input file "3d.txt"
    if ( std::string(inputFile).find("3d.txt") != std::string::npos) {
        ISPDParser::Net *net = ispdData->nets[0];
        // Decompose multi-pin nets into two-pin nets
        // Since there are only 2 pins in the given net, this step is trivial
        net->twopin.push_back(ISPDParser::TwoPin());
        ISPDParser::TwoPin &twoPin = net->twopin.back();
        twoPin.from = net->pin2D[0];
        twoPin.to   = net->pin2D[1];

        // Assume the two pin net is routed
        // The following code is to assign routing paths to the two-pin net
        // The routing path is a sequence of routing segments
        // For a horizontal segment, the start point is the left grid coordinate
        // For a vertical segment, the start point is the bottom grid coordinate
        // Please check the figures in https://www.ispd.cc/contests/08/3d.pdf
        twoPin.parNet = net;
        twoPin.path.emplace_back(0, 0, true);
        twoPin.path.emplace_back(1, 0, false);
        twoPin.path.emplace_back(0, 1, true);
        twoPin.path.emplace_back(0, 1, false);
        twoPin.path.emplace_back(0, 2, true);
        twoPin.path.emplace_back(1, 2, true);
        twoPin.path.emplace_back(2, 1, false);
        twoPin.path.emplace_back(2, 0, false);
    }

    auto graph = gr.layer_assignment();
    // Output result
    graph->output3Dresult(outputFile);

    delete graph;
    delete ispdData;

    return 0;
}
