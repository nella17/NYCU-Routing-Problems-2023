#include "router.hpp"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>
#include <random>
#include "utils.hpp"

Router::Point::Point(int _x, int _y, int _z): x(_x), y(_y), z(_z) {}

Router::Clause::Clause(): weight(0), vars{} {}
Router::Clause::Clause(std::vector<int> v): weight(0), vars(v) {}
Router::Clause::Clause(std::initializer_list<int> i): weight(0), vars{ i } {}

int& Router::Clause::emplace_back(int x) {
    assert(x != 0);
    return vars.emplace_back(x);
}

Router::Router() {}

void Router::readCircuitSpec(std::ifstream& inputFile) {
    std::string word;
    inputFile >> word >> num_z;
    inputFile >> word >> min_x >> min_y >> max_x >> max_y;
    inputFile >> word >> min_pitch_size;
    inputFile >> word >> via_cost;
    inputFile >> word >> num_nets;
    nets.resize((size_t)num_nets);
    for (int netid = 0; netid < num_nets; netid++) {
        auto &net = nets[(size_t)netid];
        net.id = netid;
        inputFile >> net.name >> net.s >> net.t;
    }

    m = (int)std::ceil(std::log2(num_nets));
}

inline void generate_coors(std::vector<int>& v, std::set<int>& s, int mn, int mx, int pitch) {
    v.clear(); v.reserve((size_t)((mx - mn) / pitch));
    v.emplace_back(mn);
    for (auto x: s) {
        while (v.back() + pitch < x)
            v.emplace_back(v.back() + pitch);
        while (!v.empty() and x - v.back() < pitch)
            v.pop_back();
        v.emplace_back(x);
    }
    while (v.back() + pitch <= mx)
        v.emplace_back(v.back() + pitch);
}

size_t Router::id(int x, int y, int z) {
    return size_t( x * num_y * num_z + y * num_z + z );
}

size_t Router::id(Point p) {
    return id(p.x, p.y, p.z);
}
    
void Router::print_node(const std::vector<int>& v) {
    for (int z = 0; z < num_z; z++) {
        for (int x = 0; x < num_x; x++) for (int y = 0; y < num_y; y++) {
            std::cerr << std::setw(std::max(2, (int)std::ceil(std::log10(num_nets))))
                << v[id(x, y, z)] << " \n"[y+1==num_y];
        }
        std::cerr _ std::string(10, '-') _ std::endl;
    }
}

void Router::generateGraph() {
    std::set<int> sX, sY{};
    for (auto &net: nets) {
        sX.emplace(net.s.x);
        sX.emplace(net.t.x);
        sY.emplace(net.s.y);
        sY.emplace(net.t.y);
    }
    generate_coors(corX, sX, min_x, max_x, min_pitch_size);
    num_x = (int)corX.size();
    generate_coors(corY, sY, min_y, max_y, min_pitch_size);
    num_y = (int)corY.size();
    num_node = num_x * num_y * num_z;
    pin_node.assign((size_t)num_node, -1);
    for (int netid = 0; netid < num_nets; netid++) {
        auto &net = nets[(size_t)netid];
        net.s.x = int( std::lower_bound(ALL(corX), net.s.x) - corX.begin() );
        net.t.x = int( std::lower_bound(ALL(corX), net.t.x) - corX.begin() );
        net.s.y = int( std::lower_bound(ALL(corY), net.s.y) - corY.begin() );
        net.t.y = int( std::lower_bound(ALL(corY), net.t.y) - corY.begin() );
        pin_node[id(net.s.x, net.s.y, net.s.z)] = netid;
        pin_node[id(net.t.x, net.t.y, net.t.z)] = netid;
    }
    print_node(pin_node);
}

std::vector<int> Router::nearedge(int x, int y, int z, int netid) {
    std::vector<int> nei{};
    nei.reserve(6);
    if (x >= 1)
        nei.emplace_back(netid + varsE[id(x-1, y, z)][DIR::X]);
    if (x+1 < num_x)
        nei.emplace_back(netid + varsE[id(x, y, z)][DIR::X]);
    if (y >= 1)
        nei.emplace_back(netid + varsE[id(x, y-1, z)][DIR::Y]);
    if (y+1 < num_y)
        nei.emplace_back(netid + varsE[id(x, y, z)][DIR::Y]);
    if (z >= 1)
        nei.emplace_back(netid + varsE[id(x, y, z-1)][DIR::Z]);
    if (z+1 < num_z)
        nei.emplace_back(netid + varsE[id(x, y, z)][DIR::Z]);
    return nei;
}

std::vector<int> Router::nearedge(Point p, int netid) {
    return nearedge(p.x, p.y, p.z, netid);
}

std::vector<Router::Point> Router::neighbor(Point p) {
    auto x = p.x, y = p.y, z = p.z;
    std::vector<Point> nei{};
    nei.reserve(6);
    if (x >= 1)
        nei.emplace_back(x-1, y, z);
    if (x+1 < num_x)
        nei.emplace_back(x+1, y, z);
    if (y >= 1)
        nei.emplace_back(x, y-1, z);
    if (y+1 < num_y)
        nei.emplace_back(x, y+1, z);
    if (z >= 1)
        nei.emplace_back(x, y, z-1);
    if (z+1 < num_z)
        nei.emplace_back(x, y, z+1);
    return nei;
}

void Router::generateClauses(std::ofstream& outputFile) {
    node0           = 1;
    edge0[DIR::X]   = node0  + num_x * num_y * num_z * (m+1);
    edge0[DIR::Y]   = edge0[DIR::X] + (num_x-1) * num_y * num_z * num_nets;
    edge0[DIR::Z]   = edge0[DIR::Y] + num_x * (num_y-1) * num_z * num_nets;
    variables       = edge0[DIR::Z] + num_x * num_y * (num_z-1) * num_nets;

    std::cerr _ node0 _ edge0[DIR::X] _ edge0[DIR::Y] _ edge0[DIR::Z] _ variables _ std::endl;

    varsN.assign((size_t)num_node, 0);
    for (int x = 0; x < num_x; x++) for (int y = 0; y < num_y; y++) for (int z = 0; z < num_z; z++)
        varsN[id(x, y, z)] = node0 + x * num_y * num_z * (m+1) + y * num_z * (m+1) + z * (m+1);

    varsE.assign((size_t)num_node, {});

    for (size_t d = 0; d < 3; d++) {
        auto isX = d==DIR::X, isY = d==DIR::Y, isZ = d==DIR::Z;
        for (int x = 0; x+isX < num_x; x++)
        for (int y = 0; y+isY < num_y; y++)
        for (int z = 0; z+isZ < num_z; z++)
            varsE[id(x, y, z)][d] = edge0[d] + x * (num_y-isY) * (num_z-isZ) * num_nets + y * (num_z-isZ) * num_nets + z * num_nets;
    }

    // X -> Y = X' or Y
    // X <-> Y = (X or Y') and (X' or Y)
    // (X and Y)' = X' or Y'

    clauses.clear();
    auto add_clause = [&](Clause c, int w = 0) {
        // std::cerr _ "  add_clause :" _ c _ std::endl;
        c.weight = w;
        clauses.emplace_back(c);
    };

    for (int x = 0; x < num_x; x++) for (int y = 0; y < num_y; y++) for (int z = 0; z < num_z; z++) {
        auto netid = pin_node[id(x, y, z)];
        auto start = varsN[id(x, y, z)];
        if (netid == -1) {
            auto ne = nearedge(x, y, z, 0);
            auto size = ne.size();
            // std::cerr _ "non-pin node" _ x _ y _ z _ ":"; for (auto k: nei) std::cerr _ k; std::cerr _ std::endl;
            // not choose 3
            for (int nid = 0; nid < num_nets; nid++) {
                for (auto&& ids: comb_id(size, 3)) {
                    Clause c;
                    for (auto i: ids) {
                        auto k = ne[i] + nid;
                        c.emplace_back(-k);
                    }
                    add_clause(c);
                }
            }
            // not choose 1
            for (int nid = 0; nid < num_nets; nid++) {
                for (size_t i = 0; i < size; i++) {
                    Clause c;
                    for (size_t j = 0; j < size; j++) {
                        auto k = ne[j] + nid;
                        if (i == j)
                            c.emplace_back(-k);
                        else
                            c.emplace_back(k);
                    }
                    add_clause(c);
                }
            }
            // used node
            Clause c;
            c.vars.reserve(1 + size * (size_t)num_nets);
            c.emplace_back(-(start + m));
            for (int nid = 0; nid < num_nets; nid++)
                for (auto n: ne)
                    c.emplace_back(n + nid);
            add_clause(c);
        } else {
            auto ne = nearedge(x, y, z, netid);
            auto size = ne.size();
            // std::cerr _ "pin node" _ x _ y _ z _ ":"; for (auto k: nei) std::cerr _ k; std::cerr _ std::endl;
            // select 1 edge
            add_clause(ne);
            // not select 2+ edge
            for (auto&& ids: comb_id(size, 2)) {
                Clause c;
                for (auto i: ids)
                    c.emplace_back(-ne[i]);
                add_clause(c);
            }
            // set node used by netid
            for (int b = 0; b < m; b++) {
                auto id = start + b;
                auto r = (netid & (1 << b)) ? 1 : -1;
                add_clause({ id * r });
            }
            add_clause({ start + m });
        }
    }

    // std::cerr _ "edge" _ std::endl;
    for (size_t d = 0; d < 3; d++) {
        auto isX = d==DIR::X, isY = d==DIR::Y, isZ = d==DIR::Z;
        for (int x = 0; x+isX < num_x; x++)
        for (int y = 0; y+isY < num_y; y++)
        for (int z = 0; z+isZ < num_z; z++) {
            auto e = varsE[id(x, y, z)][d];
            for (auto&& ids: comb_id((size_t)num_nets, 2)) {
                Clause c;
                for (auto i: ids) {
                    auto k = e + (int)i;
                    c.emplace_back(-k);
                }
                add_clause(c);
            }
        }
    }

    for (size_t d = 0; d < 3; d++) {
        auto isX = d==DIR::X, isY = d==DIR::Y, isZ = d==DIR::Z;
        for (int x = 0; x+isX < num_x; x++)
        for (int y = 0; y+isY < num_y; y++)
        for (int z = 0; z+isZ < num_z; z++) {
            auto e = varsE[id(x, y, z)][d];
            auto n1 = varsN[id(x, y, z)];
            auto n2 = varsN[id(x+isX, y+isY, z+isZ)];
            // E{1-n} -> Nm
            // = (E{1-n}' or Nm)
            // = (&(Ei') or Nm)
            for (int netid = 0; netid < num_nets; netid++) {
                auto k = e + netid;
                add_clause({ -k, n1 + m });
                add_clause({ -k, n2 + m });
                for (int b = 0; b < m; b++) {
                    auto r = (netid & (1 << b)) ? 1 : -1;
                    add_clause({ -k, r * (b + n1) });
                    add_clause({ -k, r * (b + n2) });
                }
            }
        }
    }

    weight = 1;
    for (auto &c: clauses)
        weight += c.weight;
    for (auto &c: clauses)
        if (!c.weight)
            c.weight = weight;
    
    std::mt19937 rng(0);
    std::shuffle(ALL(clauses), rng);

    std::cout << "p wcnf" _ variables-1 _ clauses.size() _ weight << std::endl;
    outputFile << "p wcnf" _ variables-1 _ clauses.size() _ weight << '\n';
    for (auto &c: clauses)
        outputFile << c << '\n';
}

std::string Router::getSysCommand(int) {
    // return "./open-wbo clause.sat > sat_result.txt; cat clause.sat sat_result.txt";
    return "./open-wbo clause.sat > sat_result.txt";
}

bool Router::readSolverResult(std::ifstream& inputFile, int) {
    std::string line;
    std::string status, vars;
    while (std::getline(inputFile, line)) {
        switch (line[0]) {
        case 's':
            status = line.substr(2);
            break;
        case 'v':
            vars = line.substr(2);
            break;
        }
    }

    std::cerr _ status _ std::endl;
    assert(!status.empty() and status != "UNSATISFIABLE");
    assignment.resize((size_t)variables);
    std::stringstream ss(vars);
    for (int x; ss >> x; )
        if (x > 0)
            assignment[(size_t)x] = true;
        else
            assignment[(size_t)-x] = false;

    node.resize((size_t)num_node);
    for (int x = 0; x < num_x; x++) for (int y = 0; y < num_y; y++) for (int z = 0; z < num_z; z++) {
        int v = 0, start = varsN[id(x, y, z)];
        for (int b = 0; b < m; b++)
            v |= assignment[size_t(start+b)] << b;
        bool used = assignment[size_t(start + m)];
        node[id(x, y, z)] = used ? v : -1;
    }
    print_node(node);

    for (int netid = 0; netid < num_nets; netid++) {
        auto& net = nets[(size_t)netid];
        auto& path = net.path;
        path.clear();
        path.emplace_back(net.s);
        std::set<Point> vis{};
        while (path.back() != net.t) {
            auto np = neighbor(path.back());
            auto ne = nearedge(path.back(), netid);
            assert(np.size() == ne.size());
            std::vector<size_t> found{};
            for (size_t i = 0, size = np.size(); i < size; i++) {
                if (node[id(np[i])] == netid and assignment[(size_t)ne[i]]) {
                    if (path.size() < 2 or path[path.size()-2] != np[i]) {
                        found.emplace_back(i);
                    }
                }
            }
            // std::cerr _ found.size() _ std::endl;
            assert(found.size() == 1);
            path.emplace_back(np[found[0]]);
        }
    }

    return true;
}

void Router::postProcess() {
}

void Router::printRoutingResult(std::ofstream& outputFile) {
    outputFile << "x_coors" _ corX.size() << '\n';
    for (auto it = corX.begin(); it != corX.end(); )
        outputFile << *it << " \n"[++it == corX.end()];
    outputFile << "y_coors" _ corY.size() << '\n';
    for (auto it = corY.begin(); it != corY.end(); )
        outputFile << *it << " \n"[++it == corY.end()];
    std::sort(ALL(nets), [&](auto a, auto b) {
        return a.id < b.id;
    });
    for (auto& net: nets) {
        outputFile << net.name _ net.path.size() << '\n';
        for (auto& p: net.path)
            outputFile << p << '\n';
    }
}

bool operator==(const Router::Point& a, const Router::Point& b) {
    return a.x == b.x and a.y == b.y and a.z == b.z;
}

bool operator!=(const Router::Point& a, const Router::Point& b) {
    return a.x != b.x or a.y != b.y or a.z != b.z;
}

std::istream& operator>>(std::istream& is, Router::Point& p) {
    return is >> p.x >> p.y >> p.z;
}

std::ostream& operator<<(std::ostream& os, const Router::Point& p) {
    return os << p.x _ p.y _ p.z;
}

std::ostream& operator<<(std::ostream& os, const Router::Clause& c) {
    os << c.weight;
    for (auto v: c.vars)
        os _ v;
    os _ 0;
    return os;
}
