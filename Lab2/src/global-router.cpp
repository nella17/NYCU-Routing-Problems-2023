#include "global-router.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>

GlobalRouter::Edge::Edge(int _cap): cap(_cap), he(1), of(0), net{}, twopins{} {}

void GlobalRouter::Edge::push(TwoPin* twopin, int min_width, int min_spacing) {
    auto [it, insert] = net.try_emplace(twopin->parNet->id, 0);
    if (insert) it->second++;
    assert(twopins.emplace(twopin).second);
    demand += std::max(twopin->parNet->minimumWidth, min_width) + min_spacing;
    of++;
}

void GlobalRouter::Edge::pop(TwoPin* twopin, int min_width, int min_spacing) {
    auto it = net.find(twopin->parNet->id);
    if (--it->second == 0)
        net.erase(it);
    assert(twopins.erase(twopin));
    demand -= std::max(twopin->parNet->minimumWidth, min_width) + min_spacing;
}

std::array<Point,4> GlobalRouter::Box::points() const {
    return {
        Point(L, U),
        Point(R, U),
        Point(L, D),
        Point(R, D),
    };
}

ld GlobalRouter::cost(const Edge& e, int k) {
    auto dah = e.he / (C[0] + C[1] * std::sqrt(k));
    auto pe = 1 + C[2] / (1 + std::exp(C[3] * (e.cap - e.demand)));
    auto be = C[4] + C[5] / std::pow(2, k);
    return (1 + dah) * pe + be;
}

ld GlobalRouter::score(const TwoPin& twopin) {
    return C[6] * twopin.overflow + C[7] * twopin.wlen();
}

int GlobalRouter::delta(const TwoPin& twopin) {
    return (int)C[8] + (int)C[9] / twopin.reroute;
}

GlobalRouter::Edge& GlobalRouter::getEdge(int x, int y, bool hori) {
    // std::cerr _ "getEdge" _ x _ y _ hori _ std::endl;
    if (hori)
        return vedges[ (size_t)x * height + (size_t)y ];
    else
        return hedges[ (size_t)x + (size_t)y * width ];
}

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData, std::array<ld, C_SIZE> _C): 
    ispdData(_ispdData), C(_C) {}

void GlobalRouter::ripup(TwoPin* twopin) {
    assert(!twopin->ripup);
    twopin->ripup = true;
    twopin->reroute++;
    for (auto rp: twopin->path)
        getEdge(rp.x, rp.y, rp.hori).pop(twopin, min_width, min_spacing);
}

void GlobalRouter::place(TwoPin* twopin) {
    assert(twopin->ripup);
    twopin->ripup = false;
    for (auto rp: twopin->path)
        getEdge(rp.x, rp.y, rp.hori).push(twopin, min_width, min_spacing);
}

Path GlobalRouter::Lshape(Point f, Point t, int k) {
    auto Lx = [&](int y, int lx, int rx, auto func) {
        if (lx > rx) std::swap(lx, rx);
        for (auto x = lx; x < rx; x++)
            func(x, y, 1);
    };
    auto Ly = [&](int x, int ly, int ry, auto func) {
        if (ly > ry) std::swap(ly, ry);
        for (auto y = ly; y < ry; y++)
            func(x, y, 0);
    };
    auto L = [&](Point p1, Point p2, auto func) {
        // std::cerr _ "L" _ p1 _ p2 _ std::endl;
        if (p1.x == p2.x) return Ly(p1.x, p1.y, p2.y, func);
        if (p1.y == p2.y) return Lx(p1.y, p1.x, p2.x, func);
    };

    if (f.y > t.y) std::swap(f, t);
    if (f.x > t.x) std::swap(f, t);

    auto Lcost = [&](Point m) {
        ld c = 0;
        auto func = [&](int x, int y, bool hori) {
            c += cost(getEdge(x, y, hori), k);
        };
        L(f, m, func);
        L(m, t, func);
        return c;
    };

    Point m1(f.x, t.y), m2(t.x, f.y);
    auto c1 = Lcost(m1), c2 = Lcost(m2);

    auto m = (c1 != c2 ? c1 < c2 : randint(2))  ? m1 : m2;

    Path path{};
    auto func = [&](int x, int y, bool hori) {
        path.emplace_back(x, y, hori);
    };
    L(f, m, func); L(m, t, func);
    return path;
}

Path GlobalRouter::HUM(Point f, Point t, int k) {
    Path path{};
    return path;
}

void GlobalRouter::route(int timeLimitSec) {
    construct_2D_grid_graph();
    net_decomposition();
    for (auto net: ispdData->nets)
        for (auto& twopin: net->twopin)
            twopins.emplace_back(&twopin);
    init_congestion();
    pattern_routing();
    for (int k = 0; k < 1; k++)
        HUM_routing(k);
    // exit(-1);
}

void GlobalRouter::construct_2D_grid_graph() {
    // construct 2D grid graph

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ALL(ispdData->nets), [&](Net *net) {

        for (auto &_pin : net->pins) {

            int x = (std::get<0>(_pin) - ispdData->lowerLeftX) / ispdData->tileWidth;
            int y = (std::get<1>(_pin) - ispdData->lowerLeftY) / ispdData->tileHeight;
            int z = std::get<2>(_pin) - 1;

            if (std::any_of(ALL(net->pin3D), [x, y, z](const auto &pin) {
                return pin.x == x && pin.y == y && pin.z == z;
            })) continue;
            net->pin3D.emplace_back(x, y, z);

            if (std::any_of(ALL(net->pin2D), [x, y](const auto &pin) { 
                return pin.x == x && pin.y == y;
            })) continue;
            net->pin2D.emplace_back(x, y);

        }

        return net->pin3D.size() > 1000 || net->pin2D.size() <= 1;

    }), ispdData->nets.end());
    ispdData->numNet = (int)ispdData->nets.size();
}

void GlobalRouter::net_decomposition() {
    for (auto net: ispdData->nets) {
        auto sz = net->pin2D.size();
        net->twopin.clear();
        net->twopin.reserve(sz - 1);
        std::vector<bool> vis(sz, false);
        std::priority_queue<std::tuple<int, size_t, size_t>> pq{};
        auto add = [&](size_t i) {
            vis[i] = true;
            auto [xi, yi, zi] = net->pin2D[i];
            for (size_t j = 0; j < sz; j++) if (!vis[j]) {
                auto [xj, yj, zj] = net->pin2D[j];
                auto d = std::abs(xi - xj) + std::abs(yi - yj);
                pq.emplace(-d, i, j);
            }
        };
        add(0);
        while (!pq.empty()) {
            auto [d, i, j] = pq.top(); pq.pop();
            if (vis[j]) continue;
            net->twopin.emplace_back(
                net->pin2D[i],
                net->pin2D[j],
                net
            );
            add(j);
        }
        assert(net->twopin.size() == sz - 1);
        assert(net->twopin.capacity() == sz - 1);
    }
}

void GlobalRouter::init_congestion() {
    width  = (size_t)ispdData->numXGrid;
    height = (size_t)ispdData->numYGrid;
    min_width = average(ispdData->minimumWidth);
    min_spacing = average(ispdData->minimumSpacing);
    auto verticalCapacity = std::accumulate(ALL(ispdData->verticalCapacity), 0);
    auto horizontalCapacity = std::accumulate(ALL(ispdData->horizontalCapacity), 0);
    vedges.assign(width * (height - 1), Edge(verticalCapacity));
    hedges.assign((width - 1) * height, Edge(horizontalCapacity));
    for (auto capacityAdj: ispdData->capacityAdjs) {
        auto [x1,y1,z1] = capacityAdj->grid1;
        auto [x2,y2,z2] = capacityAdj->grid2;
        if (z1 != z2) continue;
        auto z = (size_t)z1 - 1;
        auto lx = std::min(x1, x2), rx = std::max(x1, x2);
        auto ly = std::min(y1, y2), ry = std::max(y1, y2);
        auto dx = rx - lx, dy = ry - ly;
        if (dx + dy != 1) continue;
        auto hori = dx;
        auto& cong = getEdge(lx, ly, hori);
        auto layerCap = (hori ? ispdData->horizontalCapacity : ispdData->verticalCapacity)[z];
        cong.cap -= layerCap - capacityAdj->reducedCapacityLevel;
        // std::cerr _ dx _ dy _ "/" _ lx _ ly _ "/" _ cong.cap _ layerCap _ std::endl;
    }
    /*
    std::cerr << "horizontalCapacity\n";
    for (int j = height-1; j >= 0; j--)
        for (int i = 0; i+1 < width; i++)
            std::cerr << getCong(i, j, 1).cap << " \n"[i+2==width];
    std::cerr << "verticalCapacity\n";
    for (int j = height-2; j >= 0; j--)
        for (int i = 0; i < width; i++)
            std::cerr << getCong(i, j, 0).cap << " \n"[i+1==width];
    //*/
}

void GlobalRouter::pattern_routing() {
    for (auto twopin: twopins) {
        twopin->ripup = true;
        twopin->path = Lshape(twopin->from, twopin->to, 0);
        // std::cerr _ *twopin _ std::endl;
        place(twopin);
    }
    std::sort(twopins.begin(), twopins.end(), [&](auto a, auto b) {
        return score(*a) > score(*b);
    });
    for (auto twopin: twopins) {
        ripup(twopin);
        twopin->path = Lshape(twopin->from, twopin->to, 0);
        place(twopin);
    }
    for (auto edges: { &vedges, &hedges }) for (auto& edge: *edges)
        edge.of = 0;
    for (auto twopin: twopins)
        twopin->reroute = 0;
}

int GlobalRouter::HUM_routing(int k) {
    std::cerr _ "HUM_routing" _ k _ std::endl;
    auto start = std::chrono::steady_clock::now();

    for (auto twopin: twopins)
        twopin->overflow = 0;
    for (auto edges: { &vedges, &hedges }) for (auto& edge: *edges) {
        edge.he += edge.of;
        edge.of = 0;
        if (edge.demand > edge.cap)
            for (auto twopin: edge.twopins)
                twopin->overflow++;
    }

    int ripupcnt = 0;
    for (auto twopin: twopins)
        if (twopin->overflow)
            ripup(twopin), ripupcnt++;
    std::cerr _ "ripup" _ ripupcnt _ "twopin" _ std::endl;

    std::sort(twopins.begin(), twopins.end(), [&](auto a, auto b) {
        if (a->ripup != b->ripup)
            return a->ripup > b->ripup;
        return score(*a) > score(*b);
    });
    for (auto twopin: twopins) {
        if (twopin->ripup) {
            twopin->path = HUM(twopin->from, twopin->to, k);
            place(twopin);
        }
    }

    std::cerr.precision(2);
    std::cerr _ "HUM_routing" _ k _ std::fixed << sec_since(start) << "s" << std::endl;

    return ripupcnt;
}

LayerAssignment::Graph* GlobalRouter::layer_assignment() {
    // Assign routing layers to the two-pin net
    auto graph = new LayerAssignment::Graph;
    graph->initialLA(*ispdData, 1);
    graph->convertGRtoLA(*ispdData, true);
    graph->COLA(true);
    return graph;
}
