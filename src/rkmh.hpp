#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <math.h>
#include <algorithm>
#include <limits>
#include "murmur3.hpp"

// From Eric's https://github.com/edawson/rkmh

namespace rkmh {

// Crazy hack char table to test for canonical bases
static const int valid_dna[127] = {
        1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
        1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1
};

// Reverse complement lookup table
static char rev_arr[26] = {
        84, 66, 71, 68, 69,
        70, 67, 72, 73, 74,
        75, 76, 77, 78, 79,
        80, 81, 82, 83, 65,
        85, 86, 87, 88, 89, 90
};

typedef uint32_t hash_t;

std::vector<hash_t> hash_sequence(const char* seq,
                                  const uint64_t& len,
                                  const uint64_t& k,
                                  const uint64_t& sketch_size);

float compare(const std::vector<hash_t>& alpha, const std::vector<hash_t>& beta, const uint64_t& k);

}
