#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>

class Router {
public:
    struct Point {
        int x, y, z;
        Point(int = 0, int = 0, int = 0);
    };

    struct Net {
        std::string name;
        Point s, t;
        std::vector<Point> path;
    };

    struct Clause {
        int weight;
        std::vector<int> vars;
        Clause();
        Clause(std::vector<int>);
        int& emplace_back(int);
    };

private:
    int num_x, num_y, num_z;
    int min_x, min_y, max_x, max_y;
    int min_pitch_size;
    int via_cost;
    int num_nets;
    std::vector<Net> nets;

    int m, num_node, node0, xedge0, yedge0, zedge0, variables, weight;
    std::vector<Clause> clauses{};

    std::vector<int> pin_node, varsN, varsX, varsY, varsZ, varsE;
    std::vector<int> corX, corY, node;
    std::vector<bool> assignment;

    size_t id(int, int, int);
    size_t id(Point);
    std::vector<int> nearedge(int, int, int, int = 0);
    std::vector<int> nearedge(Point, int = 0);
    std::vector<Point> neighbor(Point);

public:
    Router();
    void readCircuitSpec(std::ifstream& inputFile);
    void generateGraph();
    void generateClauses(std::ofstream& outputFile);
    std::string getSysCommand(int);
    bool readSolverResult(std::ifstream& inputFile, int);
    void postProcess();
    void printRoutingResult(std::ofstream& outputFile);
};

bool operator==(const Router::Point&, const Router::Point&);
bool operator!=(const Router::Point&, const Router::Point&);
std::istream& operator>>(std::istream&, Router::Point&);
std::ostream& operator<<(std::ostream&, const Router::Point&);
std::ostream& operator<<(std::ostream&, const Router::Clause&);
