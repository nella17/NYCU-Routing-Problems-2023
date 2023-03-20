#include "greedy_channel_router.hpp"

#include <algorithm>
#include <iostream>

#include "utils.hpp"

GreedyChannelRouter::Node::Node(size_t _rowId, size_t _netId):
    lastColUsed(EMPTY), rowId(_rowId), netId(_netId), nextNetId(EMPTY) {}

std::ostream& operator<<(std::ostream& os, const std::vector<GreedyChannelRouter::Node>& rowEnd) {
    std::cerr _ "rowEnd" _ std::endl;
    std::cerr _ "rowId"; for (auto &n: rowEnd) std::cerr _ n.rowId; std::cerr _ std::endl;
    std::cerr _ "lastC"; for (auto &n: rowEnd) std::cerr _ n.lastColUsed; std::cerr _ std::endl;
    std::cerr _ "netId"; for (auto &n: rowEnd) std::cerr _ n.netId; std::cerr _ std::endl;
    return os;
}

GreedyChannelRouter::Path::Path(size_t _sx, size_t _sy, size_t _ex, size_t _ey):
    sx(_sx), sy(_sy), ex(_ex), ey(_ey) {}

size_t GreedyChannelRouter::Path::length() {
    return ex - sx + ey - sy;
}

GreedyChannelRouter::End::End(size_t _col): col(_col), type(0) {}

GreedyChannelRouter::Score::Score(size_t _netId):
    netId(_netId), track(0), length(0), free(0), gap(0), blocks{} {}

void GreedyChannelRouter::Score::emplace(size_t start, size_t end, size_t cnt) {
    blocks.emplace_back(start, end);
    length += end - start;
    track += cnt;
}
void GreedyChannelRouter::Score::calc(bool keep, size_t hei) {
    free = track - blocks.size() * (keep ? 1 : 0);
    if (blocks.empty())
        gap = hei;
    else
        gap = std::min(blocks.front().first, hei - blocks.back().second);
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

bool GreedyChannelRouter::operator<(const GreedyChannelRouter& r) const {
    if (height != r.height)
        return height < r.height;
    if (cost != r.cost)
        return cost < r.cost;
    if (via_used != r.via_used)
        return via_used < r.via_used;
    return ICW < r.ICW;
}

bool GreedyChannelRouter::route() {
    if (ICW == UINT_MAX or MJL == UINT_MAX or SNC == UINT_MAX)
        return false;

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
    rowEnd.reserve(ICW+2);
    height = 0;
    while (height < ICW+2)
        rowEnd.emplace_back(height++);

    size_t step = (size_t)-1;
    // TODO
    // std::cin >> step;
    // step = 6;

    for (colIdx = 0; colIdx < width; colIdx++) {
        std::cerr _ "colIdx =" _ colIdx _ height _ std::endl;

        if (step != (size_t)-1 and colIdx == step) break;

        rowEnd.front().nextNetId = netEnds[colIdx][0];
        rowEnd.back().nextNetId  = netEnds[colIdx][1];
        liveNet.clear();

        for (size_t row = 0; row < height; row++) {
            auto& n = rowEnd[row];
            if (n.netId == n.nextNetId and n.netId != Node::EMPTY) {
                netInfos[n.netId].paths.emplace_back(
                    colIdx-1, n.rowId,
                    colIdx, n.rowId
                );
            }

            n.netId = n.nextNetId;
            n.lastColUsed = Node::EMPTY;
            if (n.netId != Node::EMPTY)
                liveNet[n.netId].emplace(row);
        }

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

        std::cerr << rowEnd;
        std::cerr _ "keepNetIds";
        for (auto x: keepNetIds)
            std::cerr _ x;
        std::cerr _ std::endl;

        size_t i = 1;
        for (auto rule: rules) {
            std::cerr _ "r" << i++ _ std::endl;
            (this->*rule)();
        }

        // std::cerr << *this << std::endl;
    }

    bool extra = false;
    for (auto &n: rowEnd)
        if (n.nextNetId != Node::EMPTY) {
            extra = true;
            break;
        }

    if (extra) {
        // TODO: split-over
        return false;
    }

    std::vector<size_t> rowIdMap(height);
    for (size_t i = 0; i < height; i++) {
        rowIdMap[ rowEnd[i].rowId ] = height - 1 - i;
        // std::cerr _ rowEnd[i].rowId _ "->" _ height-1-i _ std::endl;
    }

    wire_length = via_used = 0;

    for (auto& [netId, info]: netInfos) {
        std::map<std::pair<size_t, size_t>, size_t> cnt{};
        for (auto &p: info.paths) {
            p.sy = rowIdMap[ p.sy ];
            p.ey = rowIdMap[ p.ey ];
            if (p.sy > p.ey) std::swap(p.sy, p.ey);
            // size_t nsy = rowIdMap[ p.sy ];
            // size_t ney = rowIdMap[ p.ey ];
            // p.sx = std::min(nsy, ney);
            // p.ex = std::max(nsy, ney);
            wire_length += p.length();
            if (p.sx == p.ex)
                for (size_t y = p.sy; y <= p.ey; y++)
                    cnt[{ p.sx, y }] |= 1;
            if (p.sy == p.ey)
                for (size_t x = p.sx; x <= p.ex; x++)
                    cnt[{ x, p.sy }] |= 2;
        }
        for (auto [p, c]: cnt)
            if (c == 3)
                via_used++;
    }

    cost = wire_length + via_used * via_cost;

    return true;
}

bool GreedyChannelRouter::useV(size_t netId, size_t col, size_t begin, size_t end) {
    std::cerr _ "useV" _ netId _ col _ begin _ end _ std::endl;
    if (begin == end) return false;
    if (begin > end) std::swap(begin, end);
    for (size_t row = begin; row <= end; row++)
        if (!rowEnd[row].empty(netId))
            return false;
    bool bT = rowEnd[begin].netId == netId, bF = rowEnd[begin].netId == Node::EMPTY, bb = bT | bF;
    bool eT = rowEnd[ end ].netId == netId, eF = rowEnd[ end ].netId == Node::EMPTY, ee = eT | eF;
    if (!bT and !eT) return assert(false), false;
    if ((bT and !ee) or (eT and !bb)) return false;
    // std::cerr _ bT << bF << '/' << eT << eF << std::endl;
    netInfos[netId].paths.emplace_back(
        col, rowEnd[begin].rowId,
        col, rowEnd[end].rowId
    );
    for (size_t row = begin; row <= end; row++) {
        auto& n = rowEnd[row];
        n.lastColUsed = netId;
        if (n.netId == netId) {
            n.netId = Node::EMPTY;
            liveNet[netId].erase(row);
        }
    }
    if (bT and !eT) {
        rowEnd[ end ].netId = netId;
        liveNet[netId].emplace( end );
    }
    if (!bT and eT) {
        rowEnd[begin].netId = netId;
        liveNet[netId].emplace(begin);
    }
    // std::cerr << rowEnd;
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

bool GreedyChannelRouter::placeNet(size_t netId, size_t start, size_t end, size_t target) {
    if (start > end) std::swap(start, end);
    start = std::max(start, 1ul);
    end = std::min(end, height-2);
    std::pair<long, size_t> mn{ INT_MAX, UINT_MAX };
    for (size_t row = start; row <= end; row++)
        if (rowEnd[row].netId == Node::EMPTY)
            mn = std::min(mn, std::pair{ std::abs((long)(row - target)), row });
    auto row = mn.second;
    if (row == UINT_MAX) return false;
    // std::cerr _ "placeNet" _ netId _ row _ start _ end _ target _ std::endl;
    rowEnd[row].netId = netId;
    liveNet[netId].emplace(row);
    return true;
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

    auto put = [&](size_t begin, size_t row, size_t target = UINT_MAX) {
        auto netId = rowEnd[begin].netId;
        if (netId == Node::EMPTY) {
            rowEnd[begin].lastColUsed = Node::EMPTY_USED;
        } else {
            assert(useV(netId, colIdx, begin, row));
            if (rowEnd[row].netId == netId) return;
            if (keepNet(netId) or liveNet[netId].size() >= 1) {
                if (target == UINT_MAX) target = row;
                if (!placeNet(netId, begin, row, target)) {
                    assert(false);
                    // auto m = height / 2;
                    // rowEnd.insert(rowEnd.begin() + (long)m, Node(height++));
                    // rowEnd[m].netId = netId;
                    // liveNet[netId].emplace(m);
                }
            }
        }
    };

    if (rowEnd.front().netId == rowEnd.back().netId) {
        auto netId = rowEnd.front().netId;
        put(0, height-1, keepNet(netId) ? nt2e(netInfos[netId].type) : height / 2);
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
    bool modify = false;
    do {
        modify = false;
        std::vector<Score> scores{};
        // std::cerr _ "liveNet" _ std::endl;
        // for (auto& [netId, netSet]: liveNet) {
        //     std::cerr _ netId _ ':';
        //     for (auto x: netSet) std::cerr _ x;
        //     std::cerr _ std::endl;
        // }
        for (auto& [netId, netSet]: liveNet) {
            size_t sz = netSet.size();
            if (sz < 2) continue;
            Score score(netId);
            std::vector<size_t> netRow(ALL(netSet));
            for (size_t i = 0; i < sz; i++) {
                size_t j = i+1;
                size_t row = netRow[i];
                while (j < sz and row < height-1 and rowEnd[row].empty(netId)) {
                    if (netRow[j] == row) j++;
                    row++;
                }
                if (i < --j) {
                    score.emplace(netRow[i], netRow[j], j-i+1);
                    i = j;
                }
            }
            bool keep = keepNet(netId) or score.blocks.size() > 1;
            score.calc(keep, height);
            // std::cerr _ netId _ score.track _ score.blocks.size() _ score.length _ score.gap _ std::endl;
            if (score.empty()) continue;
            scores.emplace_back(score);
        }
        std::sort(ALL(scores));
        if (!scores.empty()) {
            modify = true;
            auto& score = scores[0];
            auto netId = score.netId;
            bool kn = keepNet(netId);
            bool keep = kn or score.blocks.size() > 1;
            size_t target = kn ? nt2e(netInfos[netId].type) : (score.blocks.front().first + score.blocks.back().second) / 2;
            // std::cerr _ netId _ score.blocks.size() _ keep _ std::endl;
            for (auto [s,e]: score.blocks) {
                // std::cerr _ s <<','<< e; std::cerr _ std::endl << rowEnd;
                assert(useV(netId, colIdx, s, e));
                if (keep) assert(placeNet(netId, s, e, target));
            }
        }
    } while (modify);
}

void GreedyChannelRouter::r3() {
}

void GreedyChannelRouter::r4() {
    std::vector<std::tuple<bool, size_t, NetType, size_t>> v{};
    for (auto& [netId, netSet]: liveNet) {
        auto& info = netInfos[netId];
        if (info.ends.empty()) continue;
        v.emplace_back(
            info.type == NetType::Steady,
            info.next(),
            info.type,
            netId
        );
    }

    std::sort(ALL(v));
    for (auto [bt, col, type, netId]: v) {
        size_t end = nt2e(type);
        // std::cerr _ col _ type _ netId _ d _ end _ std::endl;
        // if (d == 0) continue;
        for (size_t i = 1; i < height-1; i++)
            if (rowEnd[i].netId == netId and rowEnd[i].empty(netId)) {
                if (i == end) continue;
                size_t d = i < end ? 1 : (size_t)-1;
                size_t j = i;
                while (j+d != end and rowEnd[j+d].empty(netId))
                    j += d;
                while (j != i and rowEnd[j].netId != Node::EMPTY)
                    j -= d;
                // std::cerr _ netId _ col _ type _ i _ j  _ '/' _ (j-i)*d _ MJL _ std::endl;
                if ((j-i)*d >= MJL) {
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
        // std::cerr << rowEnd;
        rowEnd.insert(rowEnd.begin() + (long)row, Node(height++));
        // std::cerr << rowEnd;
        std::cerr _ "insert" _ row _ height _ std::endl;
        assert(useV(netId, colIdx, begin, row));
        rowEnd[row].netId = netId;
        liveNet[netId].emplace(row);
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
    os << "height = " << r.height << std::endl
        << "ICW = " << r.ICW << std::endl
        << "MJL = " << r.MJL << std::endl
        << "SNC = " << r.SNC << std::endl
        << "wire_length = " << r.wire_length << std::endl
        << "via_used = " << r.via_used << std::endl
        << "cost = " << r.cost << std::endl;
    if (0)
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
        return os << ".V" _ p.sx _ p.sy _ p.ey;
    if (p.sy == p.ey)
        return os << ".H" _ p.sx _ p.sy _ p.ex;
    __builtin_unreachable();
}
