#pragma once

#include <vector>
#include <array>
#include <vector>
#include <iosfwd>

class GreedyChannelRouter {
public:

    struct Node {
        static const int EMPTY = 0;
        int rowId, netId;
        Node(int, int = EMPTY);
    };

    struct Path {
        int sx, sy, ex, ey;
        Path(int,int,int,int);
    };

    struct End {
        size_t col, type;
        End(size_t);
    };

    int ICW, MJL, SNC;
    size_t route();

    std::vector<std::array<int,2>> netIds;
    std::vector<std::vector<Path>> netPaths;
    std::vector<std::vector<End>> netCols;

    std::vector<std::array<bool,2>> connect;
    int colIdx;
    std::vector<Node> rowEnd;
    std::vector<bool> lastCol;

private:
    void r1(size_t);
    void r2(size_t);

    using fp = void (GreedyChannelRouter::*)(size_t);
    using RuleMap = std::vector<fp>;
    const static RuleMap rules;
};

std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter::Path& p);
