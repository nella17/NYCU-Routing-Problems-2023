#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "ispdData.h"
#include "LayerAssignment.h"
#pragma GCC diagnostic pop

#include "global-router.hpp"

#include <chrono>
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

    auto start = std::chrono::steady_clock::now();
    auto end = start + std::chrono::seconds(timeLimitSec);

    auto ispdData = parse_input(inputFile);

    GlobalRouter gr(
        ispdData,
        { 7, 4, 150, 0.3, 30, 200, 30, 1 , 5, 30 }
    );
    gr.route(end - std::chrono::seconds(120));

    auto graph = gr.layer_assignment();
    // Output result
    graph->output3Dresult(outputFile);

    delete graph;
    delete ispdData;

    return 0;
}
