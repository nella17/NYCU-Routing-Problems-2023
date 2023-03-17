#pragma once

#include <array>
#include <climits>
#include <deque>
#include <iosfwd>
#include <map>
#include <set>
#include <vector>

class GreedyChannelRouter {
public:

    struct Node {
        static const size_t EMPTY = 0;
        static const size_t EMPTY_USED = (size_t)-1;
        size_t lastColUsed;
        size_t rowId, netId, nextNetId;
        Node(size_t, size_t = EMPTY);
        inline bool empty(size_t id = EMPTY) {
            return lastColUsed == EMPTY or lastColUsed == id;
        }
    };

    struct Path {
        size_t sx, sy, ex, ey;
        Path(size_t, size_t, size_t, size_t);
    };

    struct End {
        size_t col, type;
        End(size_t);
    };

    enum NetType: size_t {
        Steady = 0,
        Raising = 1 << 0,
        Falling = 1 << 1,
    };

    struct NetInfo {
        NetType type = NetType::Steady;
        std::vector<Path> paths{};
        std::deque<End> ends{};
        inline size_t next() {
            return ends.empty() ? UINT_MAX : ends.front().col;
        }
    };

    struct Score {
        size_t netId, track, length, free, gap;
        std::vector<std::tuple<size_t, size_t, size_t>> blocks{};
        Score(size_t);
        void calc(bool, size_t);
        bool empty() const;
        bool operator<(const Score&) const;
    };

    size_t ICW, MJL, SNC;
    size_t route();

    std::vector<std::array<size_t, 2>> netEnds;

    std::map<size_t, NetInfo> netInfos;
    std::map<size_t, std::set<size_t>> liveNet; // net -> tracks
    std::set<size_t> keepNetIds;

    size_t width, height, colIdx;
    std::vector<Node> rowEnd;

private:
    bool useV(size_t, size_t, size_t, size_t);
    void fixInfo(NetInfo&);

    bool keepNet(size_t);
    bool splitNet(size_t);

    void r1();
    void r2();
    void r3();
    void r4();
    void r5();
    void r6();

    using fp = void (GreedyChannelRouter::*)();
    using RuleMap = std::vector<fp>;
    const static RuleMap rules;
};

std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter&);
std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter::Path&);
