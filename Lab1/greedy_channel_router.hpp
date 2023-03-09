#pragma once

#include <array>
#include <iosfwd>
#include <unordered_map>
#include <vector>

class GreedyChannelRouter {
public:

    struct Node {
        static const size_t EMPTY = 0;
        size_t rowId, netId;
        Node(size_t, size_t = EMPTY);
    };

    struct Path {
        size_t sx, sy, ex, ey;
        Path(size_t, size_t, size_t, size_t);
    };

    struct End {
        size_t col, type;
        End(size_t);
    };

    enum NetType {
        Steady = 0,
        Raising = 1 << 0,
        Falling = 1 << 1,
    };

    struct NetInfo: std::vector<End> {
        NetType type = NetType::Steady;
    };

    size_t ICW, MJL, SNC;
    size_t route();

    std::vector<std::array<size_t, 2>> netIds;
    std::vector<std::array<bool, 2>> connect;
    std::vector<std::vector<Path>> netPaths;

    std::unordered_map<size_t, NetInfo> netInfos;

    size_t colIdx;
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
