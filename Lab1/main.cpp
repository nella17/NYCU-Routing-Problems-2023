#include "greedy_channel_router.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>
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
        nets.erase(0);
        is.close();
    }

    size_t density = 0;
    {
        std::map<size_t, size_t> begin{}, end{};
        for (auto netId: nets) {
            begin[netId] = netEnds.size();
             end [netId] = 0;
        }
        for (size_t i = 0; i < netEnds.size(); i++) {
            for (auto x: netEnds[i]) {
                begin[x] = std::min(begin[x], i);
                 end [x] = std::max( end [x], i);
            }
        }
        std::vector<size_t> cnt(netEnds.size()+1, 0);
        for (size_t i = 1; i <= nets.size(); i++) {
            cnt[begin[i]]++;
            cnt[end[i]+1]--;
        }
        size_t cur = 0;
        std::vector<size_t> mxid{};
        for (size_t i = 0; i <= netEnds.size(); i++) {
            cur += cnt[i];
            if (cur > density) {
                density = cur;
                mxid.clear();
            }
            if (cur == density)
                mxid.emplace_back(i);
        }
        std::cout << "density = " << density << " @ ";
        for (auto i: mxid) std::cout << i << ' ';
        std::cout << std::endl;
    }

    GreedyChannelRouter router, ans;

    router.via_cost = 5;
    router.netEnds = netEnds;

    // router.ICW = 14;
    // router.MJL = 3;
    // router.SNC = 9;

    if (router.ICW != UINT_MAX) {
        auto finish = router.route();
        if (!finish) throw "route failed";
        ans = router;
    } else {
        size_t mx_SNC = std::min(20ul, netEnds.size());
        for (size_t h = density - std::min(density, 10lu); h <= ans.height; h++)
            for (size_t mjl = 1; mjl <= h / 2; mjl++)
                for (size_t snc = 1; snc <= mx_SNC; snc++)
        {
            router.ICW = h;
            router.MJL = mjl;
            // router.MJL = std::max(1lu, h / 4);
            router.SNC = snc;
            try {
                auto finish = router.route();
                if (!finish) throw "route failed";
                if (router < ans) {
                    std::cout << router << std::endl;
                    ans = router;
                }
            } catch (...) {
                // std::getchar();
                // std::cout << "h = " << h << " failed" << std::endl;
            }
        }
        if (ans.ICW == UINT_MAX) throw "route failed";
    }

    std::cout << ans;

    {
        std::ofstream os(output);
        for (auto netId: nets) {
            os << ".begin " << netId << '\n';
            for (auto p: ans.netInfos[netId].paths)
                os << p << '\n';
            os << ".end\n";
        }
        os.close();
    }

    return 0;
}
