/** $lic$
 * Copyright (C) 2012-2015 by Massachusetts Institute of Technology
 * Copyright (C) 2010-2013 by The Board of Trustees of Stanford University
 *
 * This file is part of zsim.
 *
 * zsim is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, version 2.
 *
 * If you use this software in your research, we request that you reference
 * the zsim paper ("ZSim: Fast and Accurate Microarchitectural Simulation of
 * Thousand-Core Systems", Sanchez and Kozyrakis, ISCA-40, June 2013) as the
 * source of the simulator in any publications that use this software, and that
 * you send us a citation of your work.
 *
 * zsim is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HASH_H_
#define HASH_H_

#include <stdint.h>
#include "galloc.h"

class HashFamily : public GlobAlloc {
    public:
        HashFamily() {}
        virtual ~HashFamily() {}

        virtual uint64_t hash(uint32_t id, uint64_t val) = 0;
};

class H3HashFamily : public HashFamily {
    private:
        const uint32_t numFuncs;
        uint32_t resShift;
        uint64_t* hMatrix;
    public:
        H3HashFamily(uint32_t numFunctions, uint32_t outputBits, uint64_t randSeed = 123132127);
        virtual ~H3HashFamily();
        uint64_t hash(uint32_t id, uint64_t val);
};

class SHA1HashFamily : public HashFamily {
    private:
        int numFuncs;
        int numPasses;

        //SHA1 is quite expensive and returns large blocks, so we use memoization and chunk the block to implement hash function families.
        uint64_t memoizedVal;
        uint32_t* memoizedHashes;
    public:
        explicit SHA1HashFamily(int numFunctions);
        uint64_t hash(uint32_t id, uint64_t val);
};

/* Used when we don't want hashing, just return the value */
class IdHashFamily : public HashFamily {
    public:
        inline uint64_t hash(uint32_t id, uint64_t val) {return val;}
};

#define MAX_LENGTH 64
#define subcells sbox[sbox_use]
#define subcells_inv sbox_inv[sbox_use]

typedef unsigned char          cell_t;

// Qarma64 is a lightweight tweakable block cipher.
class Qarma64HashFamily : public HashFamily {
  public:
    // Constructor.
    Qarma64HashFamily(unsigned mask);

    // Destructor.
    ~Qarma64HashFamily() {};

    // QARMA-64 encryption
    uint64_t qarma64_enc(uint64_t plaintext, uint64_t tweak, int rounds);

    // QARMA-64 decryption
    uint64_t qarma64_dec(uint64_t plaintext, uint64_t tweak, int rounds);

    uint64_t hash(uint32_t id, uint64_t val);

  protected:
    uint64_t w0 = 0x84be85ce9804e94b;
    uint64_t k0 = 0xec2802d4e0a488e9;
    unsigned setMask;

    // int sbox_use = 0;
    // Addr check_box[3] = { 0x3ee99a6c82af0c38, 0x9f5c41ec525603c9,
    //     0xbcaf6c89de930765 };
    // int sbox_use = 1;
    // Addr check_box[3] = { 0x544b0ab95bda7c3a, 0xa512dd1e4e3ec582,
    //     0xedf67ff370a483f2 };
    int sbox_use = 2;
    uint64_t check_box[3] = { 0xc003b93999b33765, 0x270a787275c48d10,
        0x5c06a7501b63b2fd };

    int m = MAX_LENGTH / 16;

    int NUM_ENC_ROUNDS = 5;

    uint64_t alpha = 0xC0AC29B7C97C50DD;
    uint64_t c[8] = { 0x0000000000000000, 0x13198A2E03707344, 0xA4093822299F31D0,
        0x082EFA98EC4E6C89, 0x452821E638D01377, 0xBE5466CF34E90C6C,
        0x3F84D5B5B5470917, 0x9216D5D98979FB1B };

    // sbox 0: lightest version, fixed points at 0, 2.
    // sbox 1: no fixed points.
    // sbox 2: lightweight sbox from prince family.
    int sbox[3][16] = {
        { 0, 14, 2, 10, 9, 15, 8, 11, 6, 4, 3, 7, 13, 12, 1, 5},
        {10, 13, 14, 6, 15, 7, 3, 5, 9, 8, 0, 12, 11, 1, 2, 4},
        {11, 6, 8, 15, 12, 0, 9, 14, 3, 7, 4, 5, 13, 2, 1, 10}};

    int sbox_inv[3][16] = {
        { 0, 14, 2, 10, 9, 15, 8, 11,
            6, 4, 3, 7, 13, 12, 1, 5},
        {10, 13, 14, 6, 15, 7, 3, 5,
            9, 8, 0, 12, 11, 1, 2, 4},
        { 5, 14, 13, 8, 10, 11, 1, 9,
            2, 6, 15, 0, 4, 12, 7, 3}};

    int t[16] = { 0, 11, 6, 13, 10, 1, 12, 7, 5, 14, 3, 8, 15, 4,  9, 2 };
    int t_inv[16] = { 0, 5, 15, 10, 13, 8, 2, 7, 11, 14, 4, 1, 6, 3, 9, 12 };
    int h[16] = { 6, 5, 14, 15, 0, 1, 2, 3, 7, 12, 13, 4, 8, 9, 10, 11 };
    int h_inv[16] = { 4, 5, 6, 7, 11, 1, 0, 8, 12, 13, 14, 15, 9, 10, 2, 3 };

    cell_t M[16] = { 0, 1, 2, 1,
        1, 0, 1, 2,
        2, 1, 0, 1,
        1, 2, 1, 0 };

    void text2cell(cell_t* cell, uint64_t is);

    uint64_t cell2text(cell_t* cell);

    uint64_t pseudo_reflect(uint64_t is, uint64_t tk);

    uint64_t forward(uint64_t is, uint64_t tk, int r);

    uint64_t backward(uint64_t is, uint64_t tk, int r);

    cell_t LFSR(cell_t x);

    cell_t LFSR_inv(cell_t x);

    uint64_t forward_update_key(uint64_t T);

    uint64_t backward_update_key(uint64_t T);
};

#endif  // HASH_H_
