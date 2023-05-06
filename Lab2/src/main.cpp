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
    auto timeLimitSec = argc < 4 ? 60 * 30 : atoi(argv[3]);

    auto start = std::chrono::steady_clock::now();
    auto end = start + std::chrono::seconds(timeLimitSec);

    auto time = std::chrono::steady_clock::now();
    auto ispdData = parse_input(inputFile);
    std::cerr << "[*] input '" << inputFile << "' costs" _ sec_since(time) << "s" << std::endl;

    GlobalRouter gr(
        ispdData,
        { 
            7, 4,
            150, 0.3,
            30, 200,
            30, 0.2 , 1,
            5, 30
        }
    );

    time = std::chrono::steady_clock::now();
    std::condition_variable cv;
    std::mutex cv_m;
    std::thread([&]() {
        std::unique_lock<std::mutex> lk(cv_m);
        auto tl = end - std::chrono::seconds(15);
        if (cv.wait_until(lk, tl) == std::cv_status::timeout)
            gr.stop = true;
    }).detach();

    try {
        gr.route();
        cv.notify_all();
    } catch (bool done) {
        if (!done)
            std::cerr _ ">>>> route stop or fail <<<<" _ std::endl;
    } catch (...) {
        std::cerr _ ">>>> unknown error <<<<" _ std::endl;
    }
    std::cerr << "[*] route costs" _ sec_since(time) << "s\n" << std::endl;

    time = std::chrono::steady_clock::now();
    auto graph = gr.layer_assignment();
    // Output result
    graph->output3Dresult(outputFile);
    std::cerr << "[*] LayerAssignment costs" _ sec_since(time) << "s\n" << std::endl;

    delete graph;
    delete ispdData;

    return 0;
}
