#include "global-router.hpp"

#pragma GCC optimize("unroll-loops")

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>
#include <chrono>
#include <iomanip>

std::ostream& operator<<(std::ostream& os, GlobalRouter::BoxCost& box) {
    for (auto y = box.U; y >= box.B; y--) {
        for (auto x = box.L; x <= box.R; x++) {
            auto p = box(x,y).from;
            if (p.has_value())
                std::cerr _ p.value();
            else
                std::cerr _ "(  S  )";
            std::cerr << std::fixed << std::setprecision(5) << box(x,y).cost;
        }
        std::cerr _ std::endl;
    }
    return os;
}

RPoint make(Point f, Point t) {
    auto dx = std::abs(f.x - t.x);
    auto dy = std::abs(f.y - t.y);
    assert(dx + dy == 1);
    if (dx == 1 and dy == 0)
        return RPoint(std::min(f.x, t.x), f.y, 1);
    if (dx == 0 and dy == 1)
        return RPoint(f.x, std::min(f.y, t.y), 0);
    __builtin_unreachable();
}

GlobalRouter::Edge::Edge(int _cap):
    cap(_cap), demand(0), he(1), of(0),
    used(0), cost(1)
{}

bool GlobalRouter::Edge::overflow() const {
    return cap < demand;
}

bool GlobalRouter::Edge::push(TwoPin* twopin) {
    if (twopin->overflow) of++;

    bool zero = used == 0;
    if (zero) demand++;
    used++;
    return zero;
}

bool GlobalRouter::Edge::pop(TwoPin* /* twopin */) {
    bool zero = used == 1;
    if (zero) demand--;
    used--;
    return zero;
}

GlobalRouter::Box::Box(Point f, Point t):
    eL(true), eR(true), eB(true), eU(true),
    L(std::min(f.x, t.x)), R(std::max(f.x, t.x)),
    B(std::min(f.y, t.y)), U(std::max(f.y, t.y)) {}

Point GlobalRouter::Box::BL() const { return Point(L, B); }
Point GlobalRouter::Box::UR() const { return Point(R, U); }
size_t GlobalRouter::Box::width()  const { return size_t(R - L + 1); }
size_t GlobalRouter::Box::height() const { return size_t(U - B + 1); }

GlobalRouter::BoxCost::BoxCost(const Box& box):
    Box(box), cost(width() * height()) {}

GlobalRouter::BoxCost::Data& GlobalRouter::BoxCost::operator()(Point p) {
    return operator()(p.x, p.y);
}
GlobalRouter::BoxCost::Data& GlobalRouter::BoxCost::operator()(int x, int y) {
    auto i = (size_t)(x - L);
    auto j = (size_t)(y - B);
    return cost.at(i * height() + j);
}

void GlobalRouter::BoxCost::trace(Path& path, Point pp) {
    auto size = path.size() + width() * height();
    while (true) {
        auto ocp = operator()(pp).from;
        if (path.size() > size) {
            std::cerr _ "path ????" _ pp _ std::endl;
            assert(false);
        }
        if (!ocp.has_value()) break;
        auto cp = ocp.value();
        path.emplace_back(make(pp, cp));
        pp = cp;
    }
}

GlobalRouter::Net::Net(ISPDParser::Net* n):
    overflow(0), overflow_twopin(0), wlen(0), reroute(0),
    score(0), cost(0), net(n), twopins{} {}

double GlobalRouter::cost(const TwoPin* twopin) const {
    double c = 0;
    for (auto rp: twopin->path)
        c += cost(rp);
    return c;
}

double GlobalRouter::cost(Point f, Point t) const {
    auto dx = std::abs(f.x - t.x);
    auto dy = std::abs(f.y - t.y);
    if (dx == 1 and dy == 0)
        return cost(std::min(f.x, t.x), f.y, 1);
    if (dx == 0 and dy == 1)
        return cost(f.x, std::min(f.y, t.y), 0);
    return INFINITY;
}

double GlobalRouter::cost(RPoint rp) const {
    return cost(getEdge(rp));
}

double GlobalRouter::cost(int x, int y, bool hori) const {
    return cost(getEdge(x, y, hori));
}

double GlobalRouter::cost(const Edge& e) const {
    return e.cost;
}

void GlobalRouter::del_cost(Net* net) {
    for (auto twopin: net->twopins)
        for (auto rp: twopin->path)
            getEdge(rp).used++;
    for (auto twopin: net->twopins)
        del_cost(twopin);
}

void GlobalRouter::del_cost(TwoPin* twopin) {
    for (auto rp: twopin->path)
        getEdge(rp).cost = 1;
}

void GlobalRouter::add_cost(Net* net) {
    for (auto twopin: net->twopins)
        for (auto rp: twopin->path)
            getEdge(rp).used--;
    for (auto twopin: net->twopins)
        add_cost(twopin);
}

void GlobalRouter::add_cost(TwoPin* twopin) {
    for (auto rp: twopin->path) {
        auto &e = getEdge(rp);
        if (e.used == 0)
            e.cost = calc_cost(e);
    }
}

double GlobalRouter::calc_cost(const Edge& e) const {
    // const Edge& e = getEdge(x, y, hori);

    // return std::exp(std::max(1, e.demand - e.cap + 1) * 2);
    auto demand = e.demand + 1;
    auto cap = e.cap;
    auto of = demand - cap;

    auto pe = get_cost_pe(of);

    if (selcost == 2) {
        auto dah = pow(e.he, 3.6) / 100;
        auto be = 200;
        auto c = (1 + dah) * pe + be;
        return c;
    }

    // if (selcost == 1)
    //     return (demand / (cap + 1.0) + pe + e.he) * (demand > cap ? 1e6 : 10);

    return pe * 10 + 200;
}

double GlobalRouter::get_cost_pe(int of) const {
    auto i = of + COSTOFF;
    if (i <= 0)
        return cost_pe[0];
    if (i >= COSTSZ)
        return cost_pe[COSTSZ-1];
    return cost_pe[i];
}

void GlobalRouter::build_cost() {
    constexpr double z = 200;
    for (int i = 0; i < COSTSZ; i++) {
        auto of = i - COSTOFF;
        switch (selcost) {
        case 0:
            cost_pe[i] = 1 + z / (1 + std::exp(-0.3 * of));
            break;
        case 1:
            cost_pe[i] = 1 + z / (1 + std::exp(-0.5 * of));
            break;
        case 2:
            cost_pe[i] = 1 + z / (1 + std::exp(-0.7 * of));
            break;
        }
        // std::cerr << cost_pe[i] << " \n"[i+1==COSTSZ];
    }
    for (auto& edge: grid)
        edge.cost = calc_cost(edge);
}

void GlobalRouter::sort_twopins() {
    std::sort(ALL(nets), [&](auto a, auto b) {
        auto sa = score(a);
        auto sb = score(b);
        return sa > sb;
    });
    for (auto net: nets)
        std::sort(ALL(net->twopins), [&](auto a, auto b) {
            auto sa = score(a);
            auto sb = score(b);
            return sa != sb ? sa < sb : a->HPWL() < b->HPWL();
        });
}

double GlobalRouter::score(const TwoPin* twopin) const {
    if (selcost == 2)
        return 60 * twopin->overflow + 1 * twopin->wlen();

    auto dx = 1 + abs(twopin->from.x - twopin->to.x);
    auto dy = 1 + abs(twopin->from.y - twopin->to.y);

    if (selcost == 1)
        return 60 * twopin->overflow + (dx * dy);

    return 100.0 / std::max(dx, dy);
}

double GlobalRouter::score(const Net* net) const {
    return 10 * net->overflow + net->overflow_twopin + 3 * log2(net->cost);
}

int GlobalRouter::delta(const TwoPin* twopin) const {
    if (twopin->reroute <= 2)
        return 5;
    if (twopin->reroute <= 6)
        return 20;
    return 15;
}

const GlobalRouter::Edge& GlobalRouter::getEdge(RPoint rp) const {
    return getEdge(rp.x, rp.y, rp.hori);
}

GlobalRouter::Edge& GlobalRouter::getEdge(RPoint rp) {
    return getEdge(rp.x, rp.y, rp.hori);
}

const GlobalRouter::Edge& GlobalRouter::getEdge(int x, int y, bool hori) const {
    return grid.at(x, y, hori);
}

GlobalRouter::Edge& GlobalRouter::getEdge(int x, int y, bool hori) {
    return grid.at(x, y, hori);
}

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData):
    stop(false), print(true), ispdData(_ispdData) {}

GlobalRouter::~GlobalRouter() {
    for (auto net: nets) delete net;
}

void GlobalRouter::ripup(TwoPin* twopin) {
    assert(!twopin->ripup);
    twopin->ripup = true;
    twopin->reroute++;
    // auto netId = twopin->parNet->id;
    // auto net = id2net.find(netId)->second;
    for (auto rp: twopin->path)
        if (getEdge(rp).pop(twopin)) {
            // net->wlen--;
        }
}

void GlobalRouter::place(TwoPin* twopin) {
    assert(twopin->ripup);
    twopin->ripup = false;
    // auto netId = twopin->parNet->id;
    // auto net = id2net.find(netId)->second;
    for (auto rp: twopin->path)
        if (getEdge(rp).push(twopin)) {
            // net->wlen++;
        }
}

void GlobalRouter::Lshape(TwoPin* twopin) {
    auto Lx = [&](int y, int L, int R, auto func) {
        if (L > R) std::swap(L, R);
        for (auto x = L; x < R; x++)
            func(x, y, 1);
    };
    auto Ly = [&](int x, int B, int U, auto func) {
        if (B > U) std::swap(B, U);
        for (auto y = B; y < U; y++)
            func(x, y, 0);
    };
    auto L = [&](Point p1, Point p2, auto func) {
        // std::cerr _ "L" _ p1 _ p2 _ std::endl;
        if (p1.x == p2.x) return Ly(p1.x, p1.y, p2.y, func);
        if (p1.y == p2.y) return Lx(p1.y, p1.x, p2.x, func);
    };

    auto& path = twopin->path;
    auto f = twopin->from, t = twopin->to;
    if (f.y > t.y) std::swap(f, t);
    if (f.x > t.x) std::swap(f, t);

    auto Lcost = [&](Point m) {
        double c = 0;
        auto func = [&](int x, int y, bool hori) {
            c += cost(x, y, hori);
        };
        L(f, m, func);
        L(m, t, func);
        return c;
    };

    Point m1(f.x, t.y), m2(t.x, f.y);
    auto c1 = Lcost(m1), c2 = Lcost(m2);

    auto m = (c1 != c2 ? c1 < c2 : randint(2)) ? m1 : m2;

    path.clear();
    auto func = [&](int x, int y, bool hori) {
        path.emplace_back(x, y, hori);
    };
    L(f, m, func); L(m, t, func);
}

void GlobalRouter::Zshape(TwoPin* twopin) {
    auto& path = twopin->path;
    auto f = twopin->from, t = twopin->to;
    if (f.y > t.y) std::swap(f, t);
    if (f.x > t.x) std::swap(f, t);

    BoxCost boxH(Box(f, t));
    boxH(f) = {
        .cost = 0,
        .from = std::nullopt,
    };
    auto boxV = boxH;

    auto dx = sign(t.x - f.x), dy = sign(t.y - f.y);

    calcX(boxH, f.y, f.x, t.x);
    for (auto px = f.x, x = px+dx; x != t.x+dx; px = x, x += dx)
        calcY(boxH, x, f.y, t.y);
    calcX(boxH, t.y, f.x, t.x);

    calcY(boxV, f.x, f.y, t.y);
    for (auto py = f.y, y = py+dy; y != t.y+dy; py = y, y += dy)
        calcX(boxV, y, f.x, t.x);
    calcY(boxV, t.x, f.y, t.y);

    auto& box = boxV(t).cost < boxH(t).cost ? boxV : boxH;

    path.clear();
    box.trace(path, t);
}

void GlobalRouter::calcX(BoxCost& box, int y, int bx, int ex) {
    auto dx = sign(ex - bx);
    if (dx == 0) return;
    auto pc = box(bx, y).cost;
    for (auto px = bx, x = px+dx; x != ex+dx; px = x, x += dx) {
        auto cc = pc + cost(std::min(x, px), y, 1);
        auto& data = box(x, y);
        if (data.cost <= cc) {
            pc = data.cost;
        } else {
            pc = cc;
            data = {
                .cost = cc,
                .from = Point(px, y),
            };
        }
    }
}

void GlobalRouter::calcY(BoxCost& box, int x, int by, int ey) {
    auto dy = sign(ey - by);
    if (dy == 0) return;
    auto pc = box(x, by).cost;
    for (auto py = by, y = py+dy; y != ey+dy; py = y, y += dy) {
        auto cc = pc + cost(x, std::min(y, py), 0);
        auto& data = box(x, y);
        if (data.cost <= cc) {
            pc = data.cost;
        } else {
            pc = cc;
            data = {
                .cost = cc,
                .from = Point(x, py),
            };
        }
    }
}

void GlobalRouter::monotonic(TwoPin* twopin) {
    auto& path = twopin->path;
    auto f = twopin->from, t = twopin->to;
    if (f.y > t.y) std::swap(f, t);
    if (f.x > t.x) std::swap(f, t);

    BoxCost box(Box(f, t));
    box(f) = {
        .cost = 0,
        .from = std::nullopt,
    };
    calcX(box, f.y, f.x, t.x);
    calcY(box, f.x, f.y, t.y);
    auto dy = sign(t.y - f.y);
    for (auto py = f.y, y = py+dy; y != t.y+dy; py = y, y += dy) {
        for (auto px = f.x, x = px+1; x <= t.x; px = x, x++) {
            auto cx = box(x, py).cost + cost(x, std::min(y, py), 0);
            auto cy = box(px, y).cost + cost(px, y, 1);
            auto sx = cx != cy ? cx < cy : randint(2);
            if (sx)
                box(x, y) = {
                    .cost = cx,
                    .from = Point(x, py),
                };
            else
                box(x, y) = {
                    .cost = cy,
                    .from = Point(px, y),
                };
        }
    }

    path.clear();
    box.trace(path, t);
}

void GlobalRouter::VMR_impl(Point f, Point t, BoxCost& box) {
    box(f) = {
        .cost = 0,
        .from = std::nullopt,
    };
    calcX(box, f.y, box.L, box.R);
    calcX(box, f.y, box.R, box.L);
    auto dy = sign(t.y - f.y);
    for (auto py = f.y, y = py+dy; y != t.y+dy; py = y, y += dy) {
        for (auto x = box.L; x <= box.R; x++)
            box(x, y) = {
                .cost = box(x, py).cost + cost(x, std::min(y, py), 0),
                .from = Point(x, py),
            };
        calcX(box, y, box.L, box.R);
        calcX(box, y, box.R, box.L);
    }
}

void GlobalRouter::HMR_impl(Point f, Point t, BoxCost& box) {
    box(f) = {
        .cost = 0,
        .from = std::nullopt,
    };
    calcY(box, f.x, box.B, box.U);
    calcY(box, f.x, box.U, box.B);
    auto dx = sign(t.x - f.x);
    for (auto px = f.x, x = px+dx; x != t.x+dx; px = x, x += dx) {
        for (auto y = box.B; y <= box.U; y++)
            box(x, y) = {
                .cost = box(px, y).cost + cost(std::min(x, px), y, 1),
                .from = Point(px, y),
            };
        calcY(box, x, box.B, box.U);
        calcY(box, x, box.U, box.B);
    }
}

void GlobalRouter::HUM(TwoPin* twopin) {
    bool insert = false;
    if (twopin->box == nullptr) {
        insert = true;
        twopin->box = new Box(twopin->from, twopin->to);
    }
    auto& box = *(Box*)twopin->box;

    if (insert or true) {
        // Congestion-aware Bounding Box Expansion
        std::array<int,2> CntOE{ 0, 0 };
        for (auto rp: twopin->path)
            if (getEdge(rp).overflow())
                CntOE[ rp.hori ] ++;
        auto d = delta(twopin);
        auto [cV, cH] = CntOE;
        auto lr = cV != cH ? cV > cH : randint(2); // TODO
        if ((lr and box.width() != width-1) or (!lr and box.height() == height-1)) {
            if (box.eL) box.L = std::max(0, box.L - d);
            if (box.eR) box.R = std::min((int)width-1, box.R + d);
        } else {
            if (box.eB) box.B = std::max(0, box.B - d);
            if (box.eU) box.U = std::min((int)height-1, box.U + d);
        }
    } else {
        auto d = delta(twopin);
        if (box.eL) box.L = std::max(0, box.L - d);
        if (box.eR) box.R = std::min((int)width-1, box.R + d);
        if (box.eB) box.B = std::max(0, box.B - d);
        if (box.eU) box.U = std::min((int)height-1, box.U + d);
    }

    auto& path = twopin->path;
    auto f = twopin->from, t = twopin->to;

    BoxCost CostVF(box), CostHF(box), CostVT(box), CostHT(box);
    if (std::abs(f.x - t.x) == box.R - box.L) {
        VMR_impl(f, box.BL(), CostVF); VMR_impl(f, box.UR(), CostVF);
        VMR_impl(t, box.BL(), CostVT); VMR_impl(t, box.UR(), CostVT);
    } else if (std::abs(f.y - t.y) == box.U - box.B) {
        HMR_impl(f, box.BL(), CostHF); HMR_impl(f, box.UR(), CostHF);
        HMR_impl(t, box.BL(), CostHT); HMR_impl(t, box.UR(), CostHT);
    } else {
        VMR_impl(f, box.BL(), CostVF); VMR_impl(f, box.UR(), CostVF);
        HMR_impl(f, box.BL(), CostHF); HMR_impl(f, box.UR(), CostHF);
        VMR_impl(t, box.BL(), CostVT); VMR_impl(t, box.UR(), CostVT);
        HMR_impl(t, box.BL(), CostHT); HMR_impl(t, box.UR(), CostHT);
    }
#ifdef DEBUG
    //*
    std::cerr _ f _ t _ " - " _ box.BL() _ box.UR() _ std::endl;
    std::cerr
        _ "CostVF\n" << CostVF
        _ "CostHF\n" << CostHF
        _ "CostVT\n" << CostVT
        _ "CostHT\n" << CostHT
        ; //*/
    print_edges(box.L, box.R, box.B, box.U);
#endif
    auto cF = [&](int x, int y) {
        return std::min(CostVF(x,y).cost, CostHF(x,y).cost);
    };
    auto cT = [&](int x, int y) {
        return std::min(CostVT(x,y).cost, CostHT(x,y).cost);
    };
    auto calc = [&](int x, int y) {
        return cF(x, y) + cT(x, y);
    };
    auto mx = box.L, my = box.B;
    auto mc = calc(mx, my);
    for (auto x = box.L; x <= box.R; x++) for (auto y = box.B; y <= box.U; y++) {
        auto c = calc(x, y);
        if (c < mc)
            mx = x, my = y, mc = c;
        // std::cerr << c << " \n"[y==box.U];
    }
    // std::cerr _ mx _ my _ std::endl;
    path.clear();
    Point m(mx, my);
    auto trace = [&](BoxCost& CostV, BoxCost& CostH) {
        auto& cost = CostV(mx, my).cost < CostH(mx, my).cost ? CostV : CostH;
        cost.trace(path, m);
    };
    trace(CostVF, CostHF);
    trace(CostVT, CostHT);

    constexpr double alpha = 1;
    auto update = [&](int L, int R, int B, int U) {
        auto ec = calc(L, B);
        for (int ux = L; ux <= R; ux++) for (int uy = B; uy <= U; uy++)
            for (int vx = L; vx <= R; vx++) for (int vy = B; vy <= U; vy++) {
                auto d = std::abs(ux - vx) + std::abs(uy - vy);
                auto c = cF(ux, uy) + cT(vx, vy) + d * alpha;
                if (c < ec) ec = c;
            }
        return mc >= ec;
    };
    box.eL = update(box.L, box.L, box.B, box.U);
    box.eR = update(box.R, box.R, box.B, box.U);
    box.eB = update(box.L, box.R, box.B, box.B);
    box.eU = update(box.L, box.R, box.U, box.U);
}

void GlobalRouter::route(bool leave) {
    width  = (size_t)ispdData->numXGrid;
    height = (size_t)ispdData->numYGrid;
    min_width = average(ispdData->minimumWidth);
    min_spacing = average(ispdData->minimumSpacing);
    min_net = min_width + min_spacing;
    construct_2D_grid_graph();
    net_decomposition();

    for (auto net: nets) delete net;
    nets.clear();
    nets.reserve(ispdData->nets.size());
    auto twopin_count = std::accumulate(ALL(ispdData->nets), 0u, [&](auto s, auto net) {
        return s + net->twopin.size();
    });
    twopins.reserve(twopin_count);
    // id2net.reserve(2lu << int(log2((double)ispdData->nets.size())));
    for (auto net: ispdData->nets) {
        auto mynet = new Net(net);
        // id2net[ net->id ] = mynet;
        mynet->twopins.reserve(net->twopin.size());
        nets.emplace_back(mynet);
        for (auto& twopin: net->twopin) {
            twopins.emplace_back(&twopin);
            mynet->twopins.emplace_back(&twopin);
        }
    }

    selcost = 0;
    preroute();
    if (leave) return;
    selcost = 0;
    routing("Lshape", &GlobalRouter::Lshape, 1);
    routing("Zshape", &GlobalRouter::Zshape, 2);
    selcost = 1;
    routing("monotonic", &GlobalRouter::monotonic, 2);
    for (auto twopin: twopins)
        twopin->reroute = 0;
    selcost = 2;
    routing("HUM", &GlobalRouter::HUM, INT_MAX);
}

void GlobalRouter::construct_2D_grid_graph() {
    // construct 2D grid graph

    // Convert XY coordinates to grid coordinates
    // Delete nets that have more than 1000 sinks
    // Delete nets that have all pins inside the same tile
    ispdData->nets.erase(std::remove_if(ALL(ispdData->nets), [&](auto net) {

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

    auto verticalCapacity = std::accumulate(ALL(ispdData->verticalCapacity), 0);
    auto horizontalCapacity = std::accumulate(ALL(ispdData->horizontalCapacity), 0);
    verticalCapacity /= min_net;
    horizontalCapacity /= min_net;
    mx_cap = std::max(verticalCapacity, horizontalCapacity);
    grid.init(width, height, Edge(verticalCapacity), Edge(horizontalCapacity));
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
        auto& edge = getEdge(lx, ly, hori);
        auto layerCap = (hori ? ispdData->horizontalCapacity : ispdData->verticalCapacity)[z];
        edge.cap -= (layerCap - capacityAdj->reducedCapacityLevel) / min_net;
        // std::cerr _ dx _ dy _ "/" _ lx _ ly _ "/" _ cong.cap _ layerCap _ std::endl;
    }
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

void GlobalRouter::print_edges(int L, int R, int B, int U) {
    if (L == -1) L = 0;
    if (R == -1) R = (int)width-1;
    if (B == -1) B = 0;
    if (U == -1) U = (int)height-1;

    std::cerr << "horizontal edges\n";
    for (int j = U; j >= B; j--) {
        for (int i = L; i+1 <= R; i++) {
            auto& e = getEdge(i, j, 1);
            std::cerr _ '(' << i << ',' << j << ')'
                << e.demand <<'/'<< e.cap << '(' << calc_cost(e) << ')';
        }
        std::cerr _ std::endl;
    }
    std::cerr << "vertical edges\n";
    for (int j = U-1; j >= B; j--) {
        for (int i = L; i <= R; i++) {
            auto& e = getEdge(i, j, 0);
            std::cerr _ '(' << i << ',' << j << ')'
                << e.demand <<'/'<< e.cap << '(' << calc_cost(e) << ')';
        }
        std::cerr _ std::endl;
    }
}

void GlobalRouter::preroute() {
    k = 0;
    if (print) std::cerr << "[*]" _ "preroute" _ std::endl;
    auto start = std::chrono::steady_clock::now();
    for (auto net: nets) {
        net->wlen = 0;
        for (auto twopin: net->twopins)
            net->wlen += twopin->HPWL();
    }
    sort_twopins();
    build_cost();
    for (auto net: nets) {
        net->wlen = 0;
        for (auto twopin: net->twopins) {
            twopin->ripup = true;
            Lshape(twopin);
            place(twopin);
            del_cost(twopin);
        }
        add_cost(net);
    }
    if (print) std::cerr _ "time" _ sec_since(start) << 's';
    check_overflow();
}

int GlobalRouter::check_overflow() {
    int mxof = 0, totof = 0;

    for (auto& edge: grid) {
        edge.he += edge.of;
        edge.of = 0;
        if (edge.overflow()) {
            auto of = edge.demand - edge.cap;
            totof += of;
            if (of > mxof) mxof = of;
        }
    }

    int ofnet = 0, oftp = 0, wl = 0;

    for (auto net: nets) {
        net->cost = net->wlen = net->overflow = net->overflow_twopin = 0;
        for (auto twopin: net->twopins) {
            twopin->overflow = 0;
            for (auto rp: twopin->path) {
                auto& e = getEdge(rp);
                bool zero = e.used++ == 0;
                if (zero) net->wlen++;
                if (e.overflow()) {
                    twopin->overflow++;
                    if (zero) {
                        net->cost += cost(e);
                        net->overflow++;
                    }
                }
            }
            if (twopin->overflow) {
                net->overflow_twopin++;
                oftp++;
            }
        }
        wl += net->wlen;
        if (net->overflow)
            ofnet++;
        for (auto twopin: net->twopins)
            for (auto rp: twopin->path)
                getEdge(rp).used--;
    }

    if (print) std::cerr 
        _ "  tot overflow" _ totof
        _ "  mx overflow" _ mxof
        _ "  wirelength" _ wl
        _ "  of net" _ ofnet
        _ "  of twopin" _ oftp
        _ std::endl;
    return totof;
}

void GlobalRouter::ripup_place(FP fp) {
    sort_twopins();
    for (auto net: nets) {
        for (auto twopin: net->twopins) {
            twopin->overflow = 0;
            for (auto rp: twopin->path)
                if (getEdge(rp).overflow()) {
                    twopin->overflow = 1;
                    break;
                }
        }

        del_cost(net);

        for (auto twopin: net->twopins) {
            if (twopin->overflow) {
                ripup(twopin);
                add_cost(twopin);
            }
        }

        for (auto twopin: net->twopins) {
            if (twopin->ripup) {
                (this->*fp)(twopin);
                place(twopin);
                del_cost(twopin);
            }
        }

        add_cost(net);
    }
    if (stop) throw false;
    return;
}

void GlobalRouter::routing(const char* name, FP fp, int iteration) {
    if (print) std::cerr << "[*]" _ name _ "routing" << std::endl;
    auto start = std::chrono::steady_clock::now();
    build_cost();
    for (int i = 1; i <= iteration; i++, k++) {
        ripup_place(fp);
        if (print) std::cerr _ i _ " time" _ sec_since(start) << 's';
        if (check_overflow() == 0) throw true;
        if (stop) throw false;
    }
    if (print) std::cerr _ name _ "routing costs" _ sec_since(start) << "s" << std::endl;
}

LayerAssignment::Graph* GlobalRouter::layer_assignment(bool print_to_screen) {
    // Assign routing layers to the two-pin net
    auto graph = new LayerAssignment::Graph;
    graph->initialLA(*ispdData, 1);
    graph->convertGRtoLA(*ispdData, print_to_screen);
    graph->COLA(print_to_screen);
    return graph;
}
