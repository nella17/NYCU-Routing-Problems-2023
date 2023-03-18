#include "greedy_channel_router.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <set>
#include <string>
#include <sstream>
#include <unistd.h>

int main(int argc, const char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s [input.txt] [output.txt]\n", argv[0]);
        return -1;
    }

#ifndef DEBUG
    close(2);
#endif

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

    router.MJL = 2;
    router.SNC = 15;

    std::pair<size_t, decltype(router.netInfos)> ans{ 18, {} };

    if (ans.first != UINT_MAX) {
        router.ICW = ans.first;
        auto height = router.route();
        ans = {
            height,
            router.netInfos
        };
    } else {
        for (size_t h = 3; h <= ans.first; h++) {
            router.ICW = h;
            // router.MJL = std::max(1lu, h / 4);
            try {
                auto height = router.route();
                std::cout << "h = " << h << " height = " << height << std::endl;
                if (height < ans.first) {
                    ans = {
                        height,
                        router.netInfos
                    };
                }
            } catch (...) {
                std::cout << "h = " << h << " failed" << std::endl;
            }
        }
    }

    {
        std::ofstream os(output);
        for (auto id: nets) {
            if (ans.second[id].paths.empty()) continue;
            os << ".begin " << id << '\n';
            for (auto p: ans.second[id].paths)
                os << p << '\n';
            os << ".end\n";
        }
        os.close();
    }

    return 0;
}
