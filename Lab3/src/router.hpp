#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>

class Router {
private:
    // You can set your data structure here, like:
    /*
    std::vector<std::vector<int>> nets_coor;
    std::vector<std::vector<int>> nets_layer;
    std::vector<std::string> nets_name;
    std::vector<std::vector<int>> nets_gridline_coor;
    std::vector<int> boundary = std::vector<int>(4);
    std::vector<int> x_coors, y_coors;
    int max_layer;
    int min_pitch_size;
    int via_cost;

    //dfs can help you find the path with the result of solvers
    //Also, it can refine the result.
    std::vector<std::vector<int>> dfs(int, std::vector<int>);
    */
public:
    Router();
    void readCircuitSpec(std::ifstream& inputFile);
    void generateGraph();
    std::string getSysCommand(int);
    void generateClauses(std::ofstream& outputFile);
    bool readSolverResult(std::ifstream& inputFile, int);
    void postProcess();
    void printRoutingResult(std::ofstream& outputFile);
};
