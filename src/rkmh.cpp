#include "rkmh.hpp"

#include <limits>
#include <algorithm>
#include <cmath>

namespace rkmh {

// Check a string (as a char*) for non-canonical DNA bases
inline bool canonical(const char *x, int len) {
    bool trip = false;
    for (int i = 0; i < len; ++i) {
        trip |= valid_dna[x[i]];
    }
    return !trip;
};

/* Reverse complement the string seq
 * (assumes seq is DNA, and returns non-ACTG letters as-is*/

/* Reverse complement a C string
 * NB: does not check safety of string lengths.
 * NB: ret is modified to hold the reverse complement of seq.
 * */

inline void reverse_complement(const char *seq, char *ret, int len) {

    //assert(seq != ret);
    if (ret == nullptr) {
        ret = new char[len + 1];
    }

    for (int i = len - 1; i >= 0; i--) {
        ret[len - 1 - i] = (char) rev_arr[(int) seq[i] - 65];
    }
    ret[len] = '\0';
};

/** Primary calc_hashes function **/
/** Takes the string to be hashed, its length,
 *  a single kmer size, a pointer to hold the hashes,
 *  and an integer to hold the number of hashes.
 *
 *  Possibly thread safe:
 *      seq, len and k are not modified
 *      new [] operator is known threadsafe
 *      User must handle hashes and numhashes properly in calling function.
 **/
inline void calc_hashes_(const char *seq, const uint64_t &len,
                         const uint64_t &k, hash_t *&hashes, int numhashes) {
    char reverse[k + 1];
    char rhash[16];
    char fhash[16];
    for (int i = 0; i < numhashes; ++i) {
        if (canonical(seq + i, k)) {
            MurmurHash3_x64_128(seq + i, k, 42, &fhash);
            hashes[i] = *((hash_t*)fhash);
            //std::cerr << "hashes[" << i << "] = " << hashes[i] << std::endl;
        } else {
            hashes[i] = std::numeric_limits<hash_t>::max();
        }
    }
};

/* Calculate all the hashes of the kmers length k of seq */
inline std::vector<hash_t> calc_hashes(const char *seq,
                                       const uint64_t& seq_length,
                                       const uint64_t& k) {
    int numhashes = seq_length - k;
    std::vector<hash_t> ret(numhashes);
    hash_t *hashes = ret.data();
    calc_hashes_(seq, seq_length, k, hashes, numhashes);
    return ret;
};

/*
inline std::vector<hash_t> calc_hashes(const char *seq, const uint64_t &len, const uint64_t& k) {
    std::vector<hash_t> ret;
    std::vector<hash_t> t = calc_hashes(seq, len, k);
    ret.insert(ret.end(), t.begin(), t.end());
    return ret;
};
*/

std::vector<hash_t> hash_sequence(const char* seq,
                                  const uint64_t& len,
                                  const uint64_t& k,
                                  const uint64_t& sketch_size) {
    std::vector<hash_t> hashes = calc_hashes(seq, len, k);
    std::sort(hashes.begin(), hashes.end());
    if (hashes.size() > sketch_size) {
        hashes.erase(hashes.begin()+sketch_size, hashes.end());
    }
    // we remove non-canonical hashes which sort last
    if (hashes.back() == std::numeric_limits<hash_t>::max()) {
        hashes.erase(std::find(hashes.begin(), hashes.end(),
                               std::numeric_limits<hash_t>::max()),
                     hashes.end());
    }
    return hashes;
}

float compare(const std::vector<hash_t>& alpha, const std::vector<hash_t>& beta, const uint64_t& k) {
    int i = 0;
    int j = 0;

    uint64_t common = 0;
    uint64_t denom = 0;

    while (i < alpha.size() && j < beta.size()) {
        if (alpha[i] == beta[j]) {
            i++;
            j++;
            common++;
        } else if (alpha[i] > beta[j]) {
            j++;
        } else {
            i++;
        }

        denom++;
    }

    // complete the union operation
    denom += alpha.size() - i;
    denom += beta.size() - j;

    float distance = 0.0;

    if (common == 0) {           // avoid inf
        distance = 1.0;
    } else if (common != denom){ // avoid -0
        //todo put a flag for denom: take the smallest between alpha.size, beta.size
        //const double jaccard = double(common) / denom;
        //distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
        //distance = -log(2 * jaccard / (1. + jaccard)) / k;
        distance = -log(2.0 * common / (double(denom) + common)) / (double)k;
        if (distance > 1) {
            distance = 1.0;
        }
    }//else {distance = 0.0;}

    return distance;
}

}
