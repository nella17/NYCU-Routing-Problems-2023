#pragma once

#include <istream>
#include <ostream>
#include <iterator>
#include <vector>
#include <string>
#include <tuple>
#include <climits>

namespace ISPDParser {

struct TwoPin;
struct RPoint;
struct Point;

using Path = std::vector<RPoint>;

class Net {

public:

    Net(const std::string &_name, int _id, int _numPins, int _minimumWidth)
              : name(_name), id(_id), numPins(_numPins), minimumWidth(_minimumWidth)  {}

    std::string name;
    int id;
    int numPins;
    int minimumWidth;
    std::vector<std::tuple<int, int, int>> pins;

    // Grid coordinates
    std::vector<Point> pin2D;
    std::vector<Point> pin3D;

    // Two pin nets
    std::vector<TwoPin> twopin;

};

class CapacityAdj {

public:
    std::tuple<int, int, int> grid1;
    std::tuple<int, int, int> grid2;
    int reducedCapacityLevel;

};

class ispdData {

public:
    ispdData() = default;
    ~ispdData() {
        for (Net *net : nets)
            delete net;
        for (CapacityAdj *capacityAdj : capacityAdjs)
            delete capacityAdj;
    }

    int numXGrid;
    int numYGrid;
    int numLayer;

    std::vector<int> verticalCapacity;
    std::vector<int> horizontalCapacity;
    std::vector<int> minimumWidth;
    std::vector<int> minimumSpacing;
    std::vector<int> viaSpacing;
    
    int lowerLeftX;
    int lowerLeftY;
    int tileWidth;
    int tileHeight;

    int numNet;
    std::vector<Net *> nets;

    int numCapacityAdj;
    std::vector<CapacityAdj *> capacityAdjs;

    friend std::ostream& operator<<(std::ostream &os, const ispdData &data) {
        os << "grid " << data.numXGrid << " " << data.numYGrid << " " << data.numLayer;
        os << "\nvertical capacity ";
        std::copy(data.verticalCapacity.begin(), data.verticalCapacity.end(), 
                std::ostream_iterator<int>(os, " "));
        os << "\nhorizontal capacity ";
        std::copy(data.horizontalCapacity.begin(), data.horizontalCapacity.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nminimum width  ";
        std::copy(data.minimumWidth.begin(), data.minimumWidth.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nminimum spacing  ";
        std::copy(data.minimumSpacing.begin(), data.minimumSpacing.end(),
                std::ostream_iterator<int>(os, " "));
        os << "\nvia spacing ";
        std::copy(data.viaSpacing.begin(), data.viaSpacing.end(),
            std::ostream_iterator<int>(os, " "));
        os << "\n" << data.lowerLeftX << " " << data.lowerLeftY << " " 
                    << data.tileWidth << " " << data.tileHeight << "\n";
        os << "num net " << data.numNet << "\n";
        for (const auto *net : data.nets) {
            os << net->name << " " << net->id << " " << net->numPins << " " << net->minimumWidth << "\n";
            for (const auto &coord : net->pins)
            os << std::get<0>(coord) << " " << std::get<1>(coord) << " " << std::get<2>(coord) << "\n";
        }
        os << data.numCapacityAdj << "\n";
        for (const auto *adj : data.capacityAdjs)
            os << std::get<0>(adj->grid1) << " " << std::get<1>(adj->grid1) << " " << std::get<2>(adj->grid1) << " "
            << std::get<0>(adj->grid2) << " " << std::get<1>(adj->grid2) << " " << std::get<2>(adj->grid2) << " "
            << adj->reducedCapacityLevel << "\n";
        return os;
    }

};

static ispdData* parse(std::istream &is) {
    ispdData *data = new ispdData();
    std::string keyword;
    is >> keyword >> data->numXGrid >> data->numYGrid >> data->numLayer;
    is >> keyword >> keyword;
    data->verticalCapacity.clear();
    int val;
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->verticalCapacity.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->horizontalCapacity.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->horizontalCapacity.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->minimumWidth.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->minimumWidth.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->minimumSpacing.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->minimumSpacing.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> keyword >> keyword;
    data->viaSpacing.clear();
    for (int i = 0; i < data->numLayer; i++) {
        is >> val;
        data->viaSpacing.push_back(val);
    }
    is.ignore(INT_MAX, '\n');
    is >> data->lowerLeftX >> data->lowerLeftY >> data->tileWidth >> data->tileHeight;
    is >> keyword >> keyword >>data->numNet;
    data->nets.clear();
    for (int i = 0; i < data->numNet; i++) {
        std::string net_name;
        int id, num_pins, min_width;
        is >> net_name >> id >> num_pins >> min_width;
        Net *net = new Net{net_name, id, num_pins, min_width};
        for (int j = 0; j < num_pins; j++) {
            int x, y, z;
            is >> x >> y >> z;
            net->pins.push_back(std::make_tuple(x, y, z));
        }
        data->nets.push_back(net);
    }
    is >> data->numCapacityAdj;
    data->capacityAdjs.clear();
    for (int i = 0; i < data->numCapacityAdj; i++) {
        int x1, y1, z1, x2, y2, z2, reduced_capacity_level;
        is >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> reduced_capacity_level;
        CapacityAdj *adj = new CapacityAdj{{x1, y1, z1}, {x2, y2, z2}, reduced_capacity_level};
        data->capacityAdjs.push_back(adj);
    }
    return data;
}


struct Point {

public:
    int x, y, z;

    Point(void) : x(0), y(0), z(0) {}
    Point(int _x, int _y) : x(_x), y(_y), z(0) {}
    Point(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}

    friend std::ostream& operator<<(std::ostream &os, const Point& p) {
        return os << '(' << p.x << ',' << p.y << ',' << p.z << ')';
    }
};

// Route Point
struct RPoint {

public:
    int x, y, z;
    bool hori;

    RPoint(): z(0) {}
    RPoint(int _x, int _y, bool _h) : x(_x), y(_y), z(0), hori(_h) {}
    RPoint(int _x, int _y, int _z, bool _h) : x(_x), y(_y), z(_z), hori(_h) {}

    friend std::ostream& operator<<(std::ostream &os, const RPoint& p) {
        return os << '(' << p.x << ',' << p.y << ',' << p.z << ',' << p.hori << ')';
    }
};

struct TwoPin
{
    Point from, to;
    Net *parNet;
    int overflow;
    bool ripup;
    int reroute;
    Path path;

    TwoPin(): TwoPin(Point(), Point()) {}

    TwoPin(Point f, Point t, Net* p = nullptr):
        from(f),
        to(t),
        parNet(p),
        overflow(0),
        ripup(false),
        reroute(0)
    {}

    int HPWL() const {
        return abs(from.x - to.x) + abs(from.y - to.y);
    }

    int wlen() const {
        return path.size() or HPWL();
    }

    friend std::ostream& operator<<(std::ostream &os, const TwoPin& twopin) {
        os << twopin.parNet->name << " " << twopin.from << " -> " << twopin.to << " : ";
        for (auto p: twopin.path)
            os << p << ' ';
        return os;
    }
};

}
