#include "greedy_channel_router.h"

#include <iostream>

GreedyChannelRouter::Node::Node(int _rowId, int _netId):
    rowId(_rowId), netId(_netId) {}

GreedyChannelRouter::Path::Path(int _sx, int _sy, int _ex, int _ey):
    sx(_sx), sy(_sy), ex(_ex), ey(_ey) {}

GreedyChannelRouter::End::End(size_t _col): col(_col), type(0) {}

size_t GreedyChannelRouter::route() {
    size_t n = netIds.size();

    netPaths.assign(n, std::vector<Path>{});
    connect.assign(n, std::array<bool,2>{ false, false });

    netCols.assign(n, std::vector<End>{});
    for (size_t c = 0; c < n; c++) {
        for (int k = 0; k < 2; k++) {
            auto netId = netIds[c][k];
            if (netId) {
                if (netCols[netId].empty() or netCols[netId].back().col != c)
                    netCols[netId].emplace_back(c);
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
    auto impl = [&](int k, int di, int begin, int end) {
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
    impl(1, -1, rowEnd.size()-1, 0);
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
