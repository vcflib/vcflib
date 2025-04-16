#include "phase.hpp"

#include "split.h"

namespace vcflib
{
    void loadPhased(std::vector<std::pair<std::string, std::string>>& haplotypes, const genotype* pop) {

        int indIndex = 0;

        for (const auto& g : pop->gts) {
            vector< string > gs = split(g, "|");
            haplotypes[indIndex].first.append(gs[0]);
            haplotypes[indIndex].first.append(gs[1]);
            indIndex += 1;
        }
    }
}