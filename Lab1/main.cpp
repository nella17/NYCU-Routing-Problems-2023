#include "greedy_channel_router.h"

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

    std::vector<std::array<int,2>> netIds{};
    std::set<int> nets;
    {
        std::ifstream is(input);
        std::string line;
        std::getline(is, line);
        std::stringstream ss(line);
        int x;
        while (ss >> x) {
            netIds.emplace_back(std::array<int,2>{ x, 0 });
            nets.emplace(x);
        }
        getline(is, line);
        ss.clear();
        ss.str(line);
        for(int i = 0; ss >> x; i++) {
            netIds[i][1] = x;
            nets.emplace(x);
        }
        is.close();
    }

    GreedyChannelRouter router;
    router.netIds = netIds;

    router.ICW = 2;
    router.MJL = 2;
    router.SNC = 3;

    auto height = router.route();
    std::cerr << "height = " << height << std::endl;

    for (auto id: nets) {
        std::cerr << "net id = " << id << std::endl;
        for (auto p: router.netPaths[id])
            std::cerr << ' ' << p << std::endl;
        std::cerr << std::endl;
    }

    {
        std::ofstream os(output);
        for (auto id: nets) {
            os << ".begin " << id << '\n';
            for (auto p: router.netPaths[id])
                os << p << '\n';
            os << ".end\n";
        }
        os.close();
    }

    return 0;
}
