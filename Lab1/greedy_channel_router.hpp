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
        size_t length();
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

    inline size_t nt2d(NetType type) {
        switch (type) {
            case NetType::Falling:
                return 1;
            case NetType::Raising:
                return (size_t)-1;
            case NetType::Steady:
                return 0;
            default:
                __builtin_unreachable();
        }
    }
    inline size_t nt2e(NetType type) {
        switch (type) {
            case NetType::Falling:
                return height - 1;
            case NetType::Raising:
                return 0;
            case NetType::Steady:
                return height / 2;
            default:
                __builtin_unreachable();
        }
    }

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
        std::vector<std::pair<size_t, size_t>> blocks{};
        Score(size_t);
        void emplace(size_t, size_t, size_t);
        void calc(bool, size_t);
        bool empty() const;
        bool operator<(const Score&) const;
    };

    size_t via_cost, wire_length, via_used, cost = UINT_MAX;
    bool operator<(const GreedyChannelRouter&) const;

    size_t ICW = UINT_MAX, MJL = UINT_MAX, SNC = UINT_MAX;
    bool route();

    std::vector<std::array<size_t, 2>> netEnds;

    std::map<size_t, NetInfo> netInfos;
    std::map<size_t, std::set<size_t>> liveNet; // net -> tracks
    std::set<size_t> keepNetIds;

    size_t width, height = UINT_MAX, colIdx;
    std::vector<Node> rowEnd;
private:
    bool useV(size_t, size_t, size_t, size_t);
    void fixInfo(NetInfo&);

    bool keepNet(size_t);
    bool placeNet(size_t, size_t, size_t, size_t);

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
