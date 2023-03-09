#include "greedy_channel_router.hpp"

#include <iostream>

GreedyChannelRouter::Node::Node(size_t _rowId, size_t _netId):
    rowId(_rowId), netId(_netId) {}

GreedyChannelRouter::Path::Path(size_t _sx, size_t _sy, size_t _ex, size_t _ey):
    sx(_sx), sy(_sy), ex(_ex), ey(_ey) {}

GreedyChannelRouter::End::End(size_t _col): col(_col), type(0) {}

size_t GreedyChannelRouter::route() {
    size_t n = netIds.size();

    connect.assign(n, std::array<bool,2>{ false, false });
    netPaths.assign(n, std::vector<Path>{});

    netInfos.clear();
    netInfos.reserve(n);
    for (size_t c = 0; c < n; c++) {
        for (size_t k = 0; k < 2; k++) {
            auto netId = netIds[c][k];
            if (netId) {
                auto& netInfo = netInfos[netId];
                if (netInfo.empty() or netInfo.back().col != c)
                    netInfo.emplace_back(c);
                netInfo.back().type |= 1 << k;
            }
        }
    }

    rowEnd.clear();
    rowEnd.reserve(ICW+2);
    for (colIdx = 0; colIdx < ICW+2; colIdx++)
        rowEnd.emplace_back(colIdx);

    for (size_t c = 0; c < n; c++) {
        lastCol.assign(rowEnd.size(), false);
        for (auto rule: rules)
            (this->*rule)(c);
    }

    size_t height = rowEnd.size();
    std::vector<size_t> rowIdMap(height);
    for (size_t i = 0; i < height; i++)
        rowIdMap[ rowEnd[i].rowId ] = height - 1 - i;
    for (auto &v: netPaths)
        for (auto &p: v) {
            p.sy = rowIdMap[ p.sy ];
            p.ey = rowIdMap[ p.ey ];
        }

    return height - 2;
}

void GreedyChannelRouter::r1(size_t c) {
    auto impl = [&](size_t k, size_t di, size_t begin, size_t end) {
        auto netId = netIds[c][k];
        if (netId == 0) {
            connect[c][k] = true;
            return;
        }

        for(auto i = begin + di; i != end; i += di) {
            auto &it = rowEnd[i];
            if (it.netId == Node::EMPTY) {
                it.netId = netId;
                connect[c][k] = true;
                netPaths[netId].emplace_back(
                    c,
                    rowEnd[begin].rowId,
                    c,
                    rowEnd[i].rowId
                );
                break;
            }
        }
    };
    impl(0,  1, 0, rowEnd.size()-1);
    impl(1, (size_t)-1, rowEnd.size()-1, 0);
}

void GreedyChannelRouter::r2(size_t c) {
}

const GreedyChannelRouter::RuleMap GreedyChannelRouter::rules{
    &GreedyChannelRouter::r1,
    &GreedyChannelRouter::r2,
};

std::ostream& operator<<(std::ostream& os, const GreedyChannelRouter::Path& p) {
    if (p.sx == p.ex)
        return os << ".V" << ' ' << p.sx << ' ' << p.sy << ' ' << p.ey;
    if (p.sy == p.ey)
        return os << ".H" << ' ' << p.sx << ' ' << p.sy << ' ' << p.ex;
    __builtin_unreachable();
}
