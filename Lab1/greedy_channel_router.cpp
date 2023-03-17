#include "greedy_channel_router.hpp"

#include <algorithm>
#include <iostream>

#include "utils.hpp"

GreedyChannelRouter::Node::Node(size_t _rowId, size_t _netId):
    lastColUsed(EMPTY), rowId(_rowId), netId(_netId), nextNetId(EMPTY) {}

GreedyChannelRouter::Path::Path(size_t _sx, size_t _sy, size_t _ex, size_t _ey):
    sx(_sx), sy(_sy), ex(_ex), ey(_ey) {}

GreedyChannelRouter::End::End(size_t _col): col(_col), type(0) {}

GreedyChannelRouter::Score::Score(size_t _netId):
    netId(_netId), track(0), length(0), free(0), gap(0), blocks{} {}

void GreedyChannelRouter::Score::calc(bool keep, size_t hei) {
    free = track - blocks.size() * (keep ? 1 : 0);
    gap = std::min(std::get<0>(blocks.front()), hei - std::get<1>(blocks.back()));
}
bool GreedyChannelRouter::Score::empty() const {
    return blocks.empty() or free == 0;
}
bool GreedyChannelRouter::Score::operator<(const Score& s) const {
    if (free != s.free)
        return free > s.free;
    if (gap != s.gap)
        return gap < s.gap;
    if (length != s.length)
        return length > s.length;
    if (blocks.size() != s.blocks.size())
        return blocks.size() > s.blocks.size();
    return netId < s.netId;
}

size_t GreedyChannelRouter::route() {
    width = netEnds.size();

    netInfos.clear();
    for (colIdx = 0; colIdx < width; colIdx++) {
        for (size_t k = 0; k < 2; k++) {
            auto netId = netEnds[colIdx][k];
            if (netId) {
                auto& netInfo = netInfos[netId];
                if (netInfo.ends.empty() or netInfo.ends.back().col != colIdx)
                    netInfo.ends.emplace_back(colIdx);
                netInfo.ends.back().type |= 1 << k;
            }
        }
    }

    keepNetIds.clear();

    for (auto& [netId,info]: netInfos)
        fixInfo(info);

    rowEnd.clear();
    rowEnd.reserve(height);
    height = 0;
    while (height < ICW+2)
        rowEnd.emplace_back(height++);

    size_t step = (size_t)-1;
    // TODO
    // std::cin >> step;
    // step = 9;

    for (colIdx = 0; colIdx < width; colIdx++) {
        std::cerr _ "colIdx =" _ colIdx _ height _ std::endl;

        if (step != (size_t)-1 and colIdx == step) break;

        rowEnd.front().nextNetId = netEnds[colIdx][0];
        rowEnd.back().nextNetId  = netEnds[colIdx][1];

        for (auto &n: rowEnd) {
            if (n.netId == n.nextNetId and n.netId != Node::EMPTY) {
                netInfos[n.netId].paths.emplace_back(
                    colIdx-1, n.rowId,
                    colIdx, n.rowId
                );
            }

            n.netId = n.nextNetId;
            n.lastColUsed = Node::EMPTY;
        }

        std::cerr _ "rowEnd.netId"; for (auto &n: rowEnd) std::cerr _ n.netId; std::cerr _ std::endl;
        std::cerr _ "rowEnd.rowId"; for (auto &n: rowEnd) std::cerr _ n.rowId; std::cerr _ std::endl;

        for (auto& [netId,info]: netInfos)
            if (info.next() == colIdx) {
                info.ends.pop_front();
                fixInfo(info);
            }

        for (auto netId: netEnds[colIdx])
            if (netId != Node::EMPTY) {
                if (netInfos[netId].ends.empty())
                    keepNetIds.erase(netId);
                else
                    keepNetIds.emplace(netId);
            }

        std::cerr _ "keepNetIds";
        for (auto x: keepNetIds)
            std::cerr _ x;
        std::cerr _ std::endl;

        size_t i = 1;
        for (auto rule: rules) {
            std::cerr _ "r" << i++ _ std::endl;
            (this->*rule)();
        }

        std::cerr << *this << std::endl;
    }

    std::vector<size_t> rowIdMap(height);
    for (size_t i = 0; i < height; i++) {
        rowIdMap[ rowEnd[i].rowId ] = height - 1 - i;
        std::cerr _ rowEnd[i].rowId _ "->" _ height-1-i _ std::endl;
    }
    for (auto& [netId, info]: netInfos)
        for (auto &p: info.paths) {
            p.sy = rowIdMap[ p.sy ];
            p.ey = rowIdMap[ p.ey ];
        }

    return height - 2;
}

bool GreedyChannelRouter::useV(size_t netId, size_t col, size_t begin, size_t end) {
    std::cerr _ "useV" _ netId _ col _ begin _ end _ std::endl;
    if (begin == end) return false;
    if (begin > end) std::swap(begin, end);
    for (size_t i = begin; i <= end; i++)
        if (!rowEnd[i].empty(netId))
            return false;
    bool bT = rowEnd[begin].netId == netId, bF = rowEnd[begin].netId == Node::EMPTY, bb = bT | bF;
    bool eT = rowEnd[ end ].netId == netId, eF = rowEnd[ end ].netId == Node::EMPTY, ee = eT | eF;
    if (!bT and !eT) return assert(false), false;
    if ((bT and !ee) or (eT and !bb)) return false;
    std::cerr _ bT << bF << '/' << eT << eF << std::endl;
    netInfos[netId].paths.emplace_back(
        col, rowEnd[begin].rowId,
        col, rowEnd[end].rowId
    );
    for (size_t i = begin; i <= end; i++) {
        auto& e = rowEnd[i];
        e.lastColUsed = netId;
        if (e.netId == netId)
            e.netId = Node::EMPTY;
    }
    std::cerr _ "rowEnd.lastC"; for (auto &n: rowEnd) std::cerr _ n.lastColUsed; std::cerr _ std::endl;
    std::cerr _ "rowEnd.netId"; for (auto &n: rowEnd) std::cerr _ n.netId; std::cerr _ std::endl;
    return true;
}

void GreedyChannelRouter::fixInfo(NetInfo& info) {
    auto end = info.ends.empty() ? 0 : info.next() + SNC;
    size_t flag = 0;
    for (auto x: info.ends)
        if (x.col < end)
            flag |= x.type;
        else
            break;
    switch (flag) {
        case NetType::Raising:
            info.type = NetType::Raising;
            break;
        case NetType::Falling:
            info.type = NetType::Falling;
            break;
        default:
            info.type = NetType::Steady;
    }
}

bool GreedyChannelRouter::keepNet(size_t netId) {
    return keepNetIds.count(netId);
}
bool GreedyChannelRouter::splitNet(size_t netId) {
    return !liveNet[netId].empty();
}

void GreedyChannelRouter::r1() {
    auto find = [&](size_t begin, size_t end, size_t di) {
        auto netId = rowEnd[begin].netId;
        if (netId == Node::EMPTY)
            return begin;
        for(auto row = begin + di; row != end; row += di) {
            auto &it = rowEnd[row];
            if (it.netId == Node::EMPTY || it.netId == netId)
                return row;
        }
        return end;
    };

    auto put = [&](size_t begin, size_t row) {
        auto netId = rowEnd[begin].netId;
        if (netId == Node::EMPTY) {
            rowEnd[begin].lastColUsed = Node::EMPTY_USED;
        } else {
            useV(netId, colIdx, begin, row);
            if (keepNet(netId))
                rowEnd[row].netId = netId;
        }
    };

    if (rowEnd.front().netId == rowEnd.back().netId) {
        auto netId = rowEnd.front().netId;
        assert(useV(netId, colIdx, 0, height-1));
    } else {
        auto mT = find(0, height-1, 1);
        auto mB = find(height-1, 0, (size_t)-1);
        if (mT < mB) {
            put(0, mT);
            put(height-1, mB);
        } else {
            if (mT-0 < height-1-mB) {
                if (mT != height-1) put(0, mT);
            } else if (mT-0 > height-1-mB) {
                if (mB != 0) put(height-1, mB);
            } else {
                // TODO
                std::cerr _ "r1" _ mT _ mB _ std::endl;
            }
        }
    }
}

void GreedyChannelRouter::r2() {
    std::map<size_t, std::vector<size_t>> netCols{};
    for(size_t row = 1; row < height-1; row++) {
        auto netId = rowEnd[row].netId ;
        if (netId) netCols[netId].emplace_back(row);
    }
    bool modify = false;
    do {
        modify = false;
        std::vector<Score> scores{};
        for (auto& [netId, netCol]: netCols) {
            size_t sz = netCol.size();
            if (sz < 2) continue;
            Score score(netId);
            for (size_t row = 1, start = 0, cnt = 0; row <= height-1; row++) {
                if (row == height-1 or !rowEnd[row].empty(netId)) {
                    if (start != 0) {
                        score.blocks.emplace_back(start, row-1, cnt);
                        score.length += row - start;
                        start = cnt = 0;
                    }
                } else {
                    if (rowEnd[row].netId == netId) {
                        cnt++;
                        score.track++;
                        if (start == 0)
                            start = row;
                    }
                }
            }
            score.calc(keepNet(netId), height);
            if (score.empty()) continue;
            scores.emplace_back(score);
            std::cerr _ netId _ score.track _ score.blocks.size() _ score.length _ score.gap _ std::endl;
        }
        std::sort(ALL(scores));
        if (!scores.empty()) {
            auto& score = scores[0];
            auto netId = score.netId;
            std::cerr _ netId _ score.blocks.size() _ std::endl;
            for (auto [s,e,c]: score.blocks) std::cerr _ s <<','<< e <<','<< c; std::cerr _ std::endl;
            for (auto [s,e,c]: score.blocks)
                useV(netId, colIdx, s, e);
        }
    } while (modify);
}

void GreedyChannelRouter::r3() {
}

void GreedyChannelRouter::r4() {
    std::vector<size_t> netIds{};
    for (auto &n: rowEnd) netIds.emplace_back(n.netId);
    std::sort(ALL(netIds));
    netIds.resize((size_t)(std::unique(ALL(netIds)) - netIds.begin()));

    std::vector<std::tuple<size_t, NetType, size_t>> v{};
    for (auto netId: netIds) {
        auto& info = netInfos[netId];
        if (info.ends.empty()) continue;
        v.emplace_back(
            info.next(),
            info.type,
            netId
        );
    }

    std::sort(ALL(v));
    for (auto [col, type, netId]: v) {
        size_t d = 0, end = 0;
        switch (type) {
            case NetType::Falling:
                d = 1;
                end = height-1;
                break;
            case NetType::Raising:
                d = (size_t)-1;
                end = 0;
                break;
            case NetType::Steady:
                continue;
            default:
                __builtin_unreachable();
        }
        for (size_t i = 1; i < height-1; i++)
            if (rowEnd[i].netId == netId and rowEnd[i].empty(netId)) {
                size_t j = i;
                while (j+d != end and rowEnd[j+d].empty(netId))
                    j += d;
                while (j != i and rowEnd[j].netId != Node::EMPTY)
                    j -= d;
                std::cerr _ netId _ col _ type _ i _ j _ std::endl;
                if (j-i >= MJL) {
                    assert(useV(netId, colIdx, i, j));
                }
            }
    }
}

void GreedyChannelRouter::r5() {
    auto impl = [&](size_t begin, size_t mid, size_t di, size_t offset = 0) {
        if (!rowEnd[begin].empty()) return;
        auto netId = rowEnd[begin].netId;
        if (netId == Node::EMPTY) return;
        size_t row = begin + di;
        while (row+di != mid and rowEnd[row+di].empty())
            row += di;
        if (begin > row) begin++;
        row += offset;
        std::cerr _ "rowEnd.rowId"; for (auto &n: rowEnd) std::cerr _ n.rowId; std::cerr _ std::endl;
        rowEnd.insert(rowEnd.begin() + (long)row, Node(height++));
        std::cerr _ "rowEnd.rowId"; for (auto &n: rowEnd) std::cerr _ n.rowId; std::cerr _ std::endl;
        std::cerr _ "rowEnd.lastC"; for (auto &n: rowEnd) std::cerr _ n.lastColUsed; std::cerr _ std::endl;
        std::cerr _ "rowEnd.netId"; for (auto &n: rowEnd) std::cerr _ n.netId; std::cerr _ std::endl;
        std::cerr _ "insert" _ row _ height _ std::endl;
        assert(useV(netId, colIdx, begin, row));
        rowEnd[row].netId = netId;
    };
    impl(0, 1 + height / 2, 1);
    impl(height-1, height / 2, (size_t)-1, 1);
}

void GreedyChannelRouter::r6() {
    // std::map<size_t, size_t> cnt{};
    for (auto &n: rowEnd) {
        n.nextNetId = n.netId;
    }
}

const GreedyChannelRouter::RuleMap GreedyChannelRouter::rules{
    &GreedyChannelRouter::r1,
    &GreedyChannelRouter::r2,
    &GreedyChannelRouter::r3,
    &GreedyChannelRouter::r4,
    &GreedyChannelRouter::r5,
    &GreedyChannelRouter::r6,
};

std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter& r) {
    for (auto& [netId, info]: r.netInfos) if (!info.paths.empty()) {
        os << "netId = " << netId << std::endl;
        for (auto p: info.paths)
            if (p.sx == p.ex)
                std::cerr << ' ' << p << std::endl;
        os << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter::Path& p) {
    if (p.sx == p.ex)
        return os << ".V" << ' ' << p.sx << ' ' << p.sy << ' ' << p.ey;
    if (p.sy == p.ey)
        return os << ".H" << ' ' << p.sx << ' ' << p.sy << ' ' << p.ex;
    __builtin_unreachable();
}
