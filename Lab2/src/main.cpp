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
#include <condition_variable>
#include <mutex>
#include <thread>

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
    auto timeLimitSec = argc < 4 ? 30 * 60 : atoi(argv[3]);

    auto start = std::chrono::steady_clock::now();
    auto end = start + std::chrono::seconds(timeLimitSec);

    auto time = std::chrono::steady_clock::now();
    auto ispdData = parse_input(inputFile);
    std::cerr << "[*] input '" << inputFile << "' costs" _ sec_since(time) << "s" << std::endl;

    GlobalRouter gr(ispdData);

    time = std::chrono::steady_clock::now();
    std::condition_variable cv;
    std::mutex cv_m;
    std::unique_lock<std::mutex> lk(cv_m);

    std::thread([&]() {
        std::chrono::duration<double> Ld = time - start;
        {
            auto TispdData = parse_input(inputFile);
            GlobalRouter Tgr(TispdData);
            Tgr.print = false;
            Tgr.route(true);
            auto Ls = std::chrono::steady_clock::now();
            auto Tgraph = Tgr.layer_assignment(false);
            Tgraph->output3Dresult(outputFile);
            Ld = std::chrono::steady_clock::now() - Ls;
            delete Tgraph;
            delete TispdData;
        }
        std::cerr << "[-] preroute LayerAssignment costs" _ Ld.count() _ "s" _ std::endl;

        auto reserve = Ld * 3 + std::chrono::seconds(5);
        auto tl = end - reserve;
        if (gr.stop or cv.wait_until(lk, tl) == std::cv_status::timeout)
            gr.stop = true;
    }).detach();

    try {
        gr.route();
    } catch (bool done) {
        if (!done)
            std::cerr _ ">>>> route stop or fail <<<<" _ std::endl;
    } catch (...) {
        std::cerr _ ">>>> unknown error <<<<" _ std::endl;
    }
    gr.stop = true;
    cv.notify_all();
    std::cerr << "[*] route costs" _ sec_since(time) << "s" << std::endl;

    time = std::chrono::steady_clock::now();
    auto graph = gr.layer_assignment();
    // Output result
    graph->output3Dresult(outputFile);
    std::cerr << "[*] LayerAssignment costs" _ sec_since(time) << "s" << std::endl;

    delete graph;
    delete ispdData;

    std::exit(EXIT_SUCCESS);

    return 0;
}
