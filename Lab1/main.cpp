#include "greedy_channel_router.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <set>
#include <string>
#include <sstream>

int main(int argc, const char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s [input.txt] [output.txt]\n", argv[0]);
        return -1;
    }

    const char* input = argv[1];
    const char* output = argv[2];

    std::vector<std::array<size_t,2>> netEnds{};
    std::set<size_t> nets;
    {
        std::ifstream is(input);
        std::string line;
        std::getline(is, line);
        std::stringstream ss(line);
        size_t x;
        while (ss >> x) {
            netEnds.emplace_back(std::array<size_t,2>{ x, 0 });
            nets.emplace(x);
        }
        getline(is, line);
        ss.clear();
        ss.str(line);
        for(size_t i = 0; ss >> x; i++) {
            netEnds[i][1] = x;
            nets.emplace(x);
        }
        is.close();
    }

    GreedyChannelRouter router;
    router.netEnds = netEnds;

    router.ICW = 2;
    router.MJL = 1;
    router.SNC = 3;

    // TODO
    router.ICW = 5;

    auto height = router.route();
    std::cerr << "height = " << height << std::endl;

    {
        std::ofstream os(output);
        for (auto id: nets) {
            os << ".begin " << id << '\n';
            for (auto p: router.netInfos[id].paths)
                os << p << '\n';
            os << ".end\n";
        }
        os.close();
    }

    return 0;
}