#include "global-router.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>
#include <chrono>

GlobalRouter::Edge::Edge(int _cap): cap(_cap), demand(0), he(1), of(0), net{}, twopins{} {}

void GlobalRouter::Edge::push(TwoPin* twopin, int minw, int mins) {
    auto [it, insert] = net.try_emplace(twopin->parNet->id, 1);
    if (!insert)
        it->second++;
    else
        demand += std::max(twopin->parNet->minimumWidth, minw) + mins;
    assert(twopins.emplace(twopin).second);
    if (twopin->overflow) of++;
}

void GlobalRouter::Edge::pop(TwoPin* twopin, int minw, int mins) {
    auto it = net.find(twopin->parNet->id);
    if (--it->second == 0) {
        net.erase(it);
        demand -= std::max(twopin->parNet->minimumWidth, minw) + mins;
    }
    assert(twopins.erase(twopin));
}

GlobalRouter::Box::Box(Point f, Point t):
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

RPoint GlobalRouter::make(Point f, Point t) {
    auto dx = std::abs(f.x - t.x);
    auto dy = std::abs(f.y - t.y);
    assert(dx + dy == 1);
    if (dx == 1 and dy == 0)
        return RPoint(std::min(f.x, t.x), f.y, 1);
    if (dx == 0 and dy == 1)
        return RPoint(f.x, std::min(f.y, t.y), 0);
    __builtin_unreachable();
}

ld GlobalRouter::cost(TwoPin* twopin)  {
    ld c = 0;
    for (auto rp: twopin->path)
        c += cost(rp);
    return c;
}

ld GlobalRouter::cost(Point f, Point t) {
    auto dx = std::abs(f.x - t.x);
    auto dy = std::abs(f.y - t.y);
    if (dx == 1 and dy == 0)
        return cost(std::min(f.x, t.x), f.y, 1);
    if (dx == 0 and dy == 1)
        return cost(f.x, std::min(f.y, t.y), 0);
    return INFINITY;
}

ld GlobalRouter::cost(RPoint rp) {
    return cost(getEdge(rp));
}

ld GlobalRouter::cost(int x, int y, bool hori) {
    return cost(getEdge(x, y, hori));
}

ld GlobalRouter::cost(const Edge& e) const {
    // return std::exp(std::max(0, e.demand - e.cap + 1) * 2);
    auto dah = pow(e.he, 2) / (C[0] + C[1] * std::sqrt(k));
    auto pe = 1 + C[2] / (1 + std::exp(C[3] * (e.cap - e.demand)));
    auto be = C[4] + C[5] / std::pow(2, k);
    return (1 + dah) * pe + be;
}

ld GlobalRouter::score(const TwoPin& twopin) const {
    return C[6] * twopin.overflow + C[7] * twopin.cost + C[8] * twopin.reroute;
}

int GlobalRouter::delta(const TwoPin& twopin) const {
    return (int)C[9] + (int)C[10] / twopin.reroute;
}

GlobalRouter::Edge& GlobalRouter::getEdge(RPoint rp) {
    return getEdge(rp.x, rp.y, rp.hori);
}

GlobalRouter::Edge& GlobalRouter::getEdge(int x, int y, bool hori) {
    // std::cerr _ "getEdge" _ x _ y _ hori _ std::endl;
    if (hori)
        return hedges.at( (size_t)x * height + (size_t)y );
    else
        return vedges.at( (size_t)x + (size_t)y * width );
}

GlobalRouter::GlobalRouter(ISPDParser::ispdData* _ispdData, std::array<ld, C_SIZE> _C): 
    stop(false), print(true), ispdData(_ispdData), C(_C) {}

void GlobalRouter::ripup(TwoPin* twopin) {
    assert(!twopin->ripup);
    twopin->ripup = true;
    twopin->reroute++;
    for (auto rp: twopin->path)
        getEdge(rp).pop(twopin, min_width, min_spacing);
}

void GlobalRouter::place(TwoPin* twopin) {
    assert(twopin->ripup);
    twopin->ripup = false;
    for (auto rp: twopin->path)
        getEdge(rp).push(twopin, min_width, min_spacing);
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
        ld c = 0;
        auto func = [&](int x, int y, bool hori) {
            c += cost(x, y, hori);
        };
        L(f, m, func);
        L(m, t, func);
        return c;
    };

    Point m1(f.x, t.y), m2(t.x, f.y);
    auto c1 = Lcost(m1), c2 = Lcost(m2);

    auto m = (c1 != c2 ? c1 < c2 : randint(2))  ? m1 : m2;

    path.clear();
    auto func = [&](int x, int y, bool hori) {
        path.emplace_back(x, y, hori);
    };
    L(f, m, func); L(m, t, func);
}

void GlobalRouter::VMR_impl(Point s, Point t, BoxCost& box) {
    auto calc = [&](int y, int bx, int ex) {
        auto dx = sign(ex - bx);
        if (dx == 0) return;
        auto pc = box(bx, y).cost;
        for (auto px = bx, x = px+dx; x != ex+dx; px = x, x += dx) {
            auto cc = pc + cost(std::min(x, px), y, 1);
            auto& data = box(x, y);
            if (data.cost <= cc) {
                cc = data.cost;
            } else {
                data = {
                    .cost = cc,
                    .from = Point(px, y),
                };
            }
            pc = cc;
        }
    };
    box(s) = {
        .cost = 0,
        .from = std::nullopt,
    };
    calc(s.y, box.L, box.R);
    calc(s.y, box.R, box.L);
    auto dy = sign(t.y - s.y);
    for (auto py = s.y, y = py+dy; y != t.y+dy; py = y, y += dy) {
        for (auto x = box.L; x <= box.R; x++)
            box(x, y) = {
                .cost = box(x, py).cost + cost(x, std::min(y, py), 0),
                .from = Point(x, py),
            };
        calc(y, box.L, box.R);
        calc(y, box.R, box.L);
    }
}

void GlobalRouter::HMR_impl(Point s, Point t, BoxCost& box) {
    auto calc = [&](int x, int by, int ey) {
        auto dy = sign(ey - by);
        if (dy == 0) return;
        auto pc = box(x, by).cost;
        for (auto py = by, y = py+dy; y != ey+dy; py = y, y += dy) {
            auto cc = pc + cost(x, std::min(y, py), 0);
            auto& data = box(x, y);
            if (data.cost <= cc) {
                cc = data.cost;
            } else {
                data = {
                    .cost = cc,
                    .from = Point(x, py),
                };
            }
            pc = cc;
        }
    };
    box(s) = {
        .cost = 0,
        .from = std::nullopt,
    };
    calc(s.x, box.B, box.U);
    calc(s.x, box.U, box.B);
    auto dx = sign(t.x - s.x);
    for (auto px = s.x, x = px+dx; x != t.x+dx; px = x, x += dx) {
        for (auto y = box.B; y <= box.U; y++)
            box(x, y) = {
                .cost = box(px, y).cost + cost(std::min(x, px), y, 1),
                .from = Point(px, y),
            };
        calc(x, box.B, box.U);
        calc(x, box.U, box.B);
    }
}

void GlobalRouter::Zshape(TwoPin* twopin) {
    // TODO
}

void GlobalRouter::monotonic(TwoPin* twopin) {
    // TODO
}

std::ostream& operator<<(std::ostream& os, GlobalRouter::BoxCost& box) {
    for (auto y = box.U; y >= box.B; y--) {
        for (auto x = box.L; x <= box.R; x++) {
            auto p = box(x,y).from;
            if (p.has_value())
                std::cerr _ p.value() << box(x,y).cost;
            else
                std::cerr _ "(  S  )       ";
        }
        std::cerr _ std::endl;
    }
    return os;
}

void GlobalRouter::HUM(TwoPin* twopin) {
    auto [it,insert] = boxs.try_emplace(twopin, twopin->from, twopin->to);
    auto& box = it->second;
    if (!insert) {
        // Congestion-aware Bounding Box Expansion
        std::array<int,2> CntOE{ 0, 0 };
        for (auto rp: twopin->path)
            CntOE[ rp.hori ] ++;
        auto d = delta(*twopin);
        auto [cV, cH] = CntOE;
        if ((cV != cH ? cV < cH : randint(2))) {
            box.L = std::max(0, box.L - d);
            box.B = std::max(0, box.B - d);
        } else {
            box.R = std::min((int)width-1, box.R + d);
            box.U = std::min((int)height-1, box.U + d);
        }
    }

    auto& path = twopin->path;
    auto f = twopin->from, t = twopin->to;

    BoxCost CostVF(box), CostHF(box), CostVT(box), CostHT(box);
    VMR_impl(f, box.BL(), CostVF); VMR_impl(f, box.UR(), CostVF);
    HMR_impl(f, box.BL(), CostHF); HMR_impl(f, box.UR(), CostHF);
    VMR_impl(t, box.BL(), CostVT); VMR_impl(t, box.UR(), CostVT);
    HMR_impl(t, box.BL(), CostHT); HMR_impl(t, box.UR(), CostHT);
#ifdef DEBUG
    //*
    std::cerr
        _ "CostVF\n" << CostVF
        _ "CostHF\n" << CostHF
        _ "CostVT\n" << CostVT
        _ "CostHT\n" << CostHT
        ; //*/
#endif
    auto calc = [&](int x, int y) {
        auto cs = std::min(CostVF(x,y).cost, CostHF(x,y).cost);
        auto ct = std::min(CostVT(x,y).cost, CostHT(x,y).cost);
        return cs + ct;
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
    auto trace = [&](BoxCost& CostV, BoxCost& CostH) {
        auto& cost = CostV(mx, my).cost < CostH(mx, my).cost ? CostV : CostH;
        Point pp(mx, my);
        while (true) {
            auto ocp = cost(pp).from;
            if (path.size() > width * height) {
                std::cerr _ "path ????" _ std::endl;
                throw false;
            }
            if (!ocp.has_value()) break;
            auto cp = ocp.value();
            path.emplace_back(make(pp, cp));
            pp = cp;
        }
    };
    trace(CostVF, CostHF);
    trace(CostVT, CostHT);
}

void GlobalRouter::route(bool preroute) {
    width  = (size_t)ispdData->numXGrid;
    height = (size_t)ispdData->numYGrid;
    min_width = average(ispdData->minimumWidth);
    min_spacing = average(ispdData->minimumSpacing);
    construct_2D_grid_graph();
    net_decomposition();
    for (auto net: ispdData->nets)
        for (auto& twopin: net->twopin)
            twopins.emplace_back(&twopin);
    init();
    if (preroute) return;
    routing("Lshape", &GlobalRouter::Lshape);
    // TODO: Zshape
    // TODO: monotonic
    routing("HUM", &GlobalRouter::HUM, INT_MAX);
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

void GlobalRouter::print_edges() {
    std::cerr << "horizontalCapacity\n";
    for (int j = (int)height-1; j >= 0; j--) {
        for (int i = 0; i+1 < (int)width; i++) {
            auto& e = getEdge(i, j, 1);
            std::cerr _ e.demand <<'/'<< e.cap << '(' << cost(e) << ')';
        }
        std::cerr _ std::endl;
    }
    std::cerr << "verticalCapacity\n";
    for (int j = (int)height-2; j >= 0; j--) {
        for (int i = 0; i < (int)width; i++) {
            auto& e = getEdge(i, j, 0);
            std::cerr _ e.demand <<'/'<< e.cap << '(' << cost(e) << ')';
        }
        std::cerr _ std::endl;
    }
}

void GlobalRouter::init() {
    if (print) std::cerr << "[*]" _ "init" _ std::endl;
    auto start = std::chrono::steady_clock::now();
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
        auto& edge = getEdge(lx, ly, hori);
        auto layerCap = (hori ? ispdData->horizontalCapacity : ispdData->verticalCapacity)[z];
        edge.cap -= layerCap - capacityAdj->reducedCapacityLevel;
        // std::cerr _ dx _ dy _ "/" _ lx _ ly _ "/" _ cong.cap _ layerCap _ std::endl;
    }
    k = 0;
    for (auto twopin: twopins) {
        twopin->ripup = true;
        Lshape(twopin);
        // std::cerr _ *twopin _ std::endl;
    }
    for (auto twopin: twopins)
        place(twopin);
    for (auto twopin: twopins)
        twopin->cost = cost(twopin);
    std::sort(ALL(twopins), [&](auto a, auto b) {
        return a->cost < b->cost;
    });
    // std::shuffle(ALL(twopins), rng);
    ripup_place(&GlobalRouter::Lshape);
    for (auto edges: { &vedges, &hedges }) for (auto& edge: *edges)
        edge.of = 0;
    for (auto twopin: twopins)
        twopin->reroute = 0;
    if (print) std::cerr _ "time" _ sec_since(start) << 's';
    check_overflow();
}

int GlobalRouter::check_overflow() {
    for (auto twopin: twopins)
        twopin->overflow = 0;
    int mxof = 0, totof = 0;
    for (auto edges: { &vedges, &hedges }) for (auto& edge: *edges) {
        edge.he += edge.of;
        edge.of = 0;
        edge.overflow = edge.demand > edge.cap;
        if (edge.overflow) {
            auto of = edge.demand - edge.cap;
            totof += of;
            if (of > mxof) mxof = of;
            for (auto twopin: edge.twopins)
                twopin->overflow++;
        }
    }

    int cnt = 0, wl = 0;
    for (auto twopin: twopins) {
        wl += twopin->wlen();
        if (twopin->overflow)
            cnt++;
    }

    if (print) std::cerr 
        _ "  tot overflow" _ totof
        _ "  mx overflow" _ mxof
        _ "  wirelength" _ wl
        _ "  of twopin" _ cnt
        _ std::endl;
    return totof;
}

void GlobalRouter::ripup_place(FP fp) {
    for (auto twopin: twopins) {
        if (stop) throw false;
        if (twopin->overflow) {
            ripup(twopin);
            (this->*fp)(twopin);
            place(twopin);
        }
    }
}

void GlobalRouter::routing(const char* name, FP fp, int iteration) {
    if (print) std::cerr << "[*]" _ name _ "routing" << std::endl;
    auto start = std::chrono::steady_clock::now();
    for (int i = 1; i <= iteration; i++, k++) {
        auto iterstart = std::chrono::steady_clock::now();
        for (auto twopin: twopins)
            twopin->cost = cost(twopin);
        std::sort(ALL(twopins), [&](auto a, auto b) {
            return score(*a) > score(*b);
        });
        // std::shuffle(ALL(twopins), rng);
        ripup_place(fp);
        if (print) std::cerr _ i _ " time" _ sec_since(iterstart) << 's';
        if (check_overflow() == 0) throw true;
#ifdef DEBUG
        print_edges();
#endif
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
