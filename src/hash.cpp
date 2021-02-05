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

#include "hash.h"
#include <stdio.h>
#include <stdlib.h>
#include "log.h"
#include "mtrand.h"

H3HashFamily::H3HashFamily(uint32_t numFunctions, uint32_t outputBits, uint64_t randSeed) : numFuncs(numFunctions) {
    MTRand rnd(randSeed);

    if (outputBits <= 8) {
        resShift = 3;
    } else if (outputBits <= 16) {
        resShift = 2;
    } else if (outputBits <= 32) {
        resShift = 1;
    } else if (outputBits <= 64) {
        resShift = 0;
    } else {
        panic("Hash function can't produce more than 64 bits of output!!");
    }

    uint32_t words = 64 >> resShift;
    hMatrix = gm_calloc<uint64_t>(words*numFuncs);
    for (uint32_t ii = 0; ii < numFuncs; ii++) {
        for (uint32_t jj = 0; jj < words; jj++) {
            uint64_t val = 0;
            for (int kk = 0; kk < 64; kk++) {
                val = val << 1;
                if (rnd.randInt() % 2 == 0) val++;
            }
            //Indeed, they are distributed around 32, but we might get better mileage by forcing 32b...
            //info("H3: Function %d Matrix 64-bit word %d has %d 1s", ii, jj, __builtin_popcountll(val));
            //if (__builtin_popcountll(val) != 32) {jj--; continue;} // no difference
            hMatrix[ii*words + jj] = val;
        }
    }
}

H3HashFamily::~H3HashFamily() {
    gm_free(hMatrix);
}

/* NOTE: This is fairly well hand-optimized. Go to the commit logs to see the speedup of this function. Main things:
 * 1. resShift indicates how many bits of output are computed (64, 32, 16, or 8). With less than 64 bits, several rounds are folded at the end.
 * 2. The output folding does not mask, the output is expected to be masked by caller.
 * 3. The main loop is hand-unrolled and optimized for ILP.
 * 4. Pre-computing shifted versions of the input does not help, as it increases register pressure.
 *
 * For reference, here is the original, simpler code (computes a 64-bit hash):
 * for (uint32_t x = 0; x < 64; x++) {
 *     res ^= val & hMatrix[id*64 + x];
 *     res = (res << 1) | (res >> 63);
 * }
 */
uint64_t H3HashFamily::hash(uint32_t id, uint64_t val) {
    uint64_t res = 0;
    assert(id >= 0 && id < numFuncs);

    // 8-way unrolled loop
    uint32_t maxBits = 64 >> resShift;
    for (uint32_t x = 0; x < maxBits; x+=8) {
        uint32_t base = (id << (6 - resShift)) + x;
        uint64_t res0 = val & hMatrix[base];
        uint64_t res1 = val & hMatrix[base+1];
        uint64_t res2 = val & hMatrix[base+2];
        uint64_t res3 = val & hMatrix[base+3];

        uint64_t res4 = val & hMatrix[base+4];
        uint64_t res5 = val & hMatrix[base+5];
        uint64_t res6 = val & hMatrix[base+6];
        uint64_t res7 = val & hMatrix[base+7];

        res ^= res0 ^ ((res1 << 1) | (res1 >> 63)) ^ ((res2 << 2) | (res2 >> 62)) ^ ((res3 << 3) | (res3 >> 61));
        res ^= ((res4 << 4) | (res4 >> 60)) ^ ((res5 << 5) | (res5 >> 59)) ^ ((res6 << 6) | (res6 >> 58)) ^ ((res7 << 7) | (res7 >> 57));
        res = (res << 8) | (res >> 56);
    }

    // Fold bits to match output
    switch (resShift) {
        case 0: //64-bit output
            break;
        case 1: //32-bit output
            res = (res >> 32) ^ res;
            break;
        case 2: //16-bit output
            res = (res >> 32) ^ res;
            res = (res >> 16) ^ res;
            break;
        case 3: //8-bit output
            res = (res >> 32) ^ res;
            res = (res >> 16) ^ res;
            res = (res >> 8) ^ res;
            break;
    }

    //info("0x%lx", res);

    return res;
}

#if _WITH_POLARSSL_

#include "polarssl/sha1.h"

SHA1HashFamily::SHA1HashFamily(int numFunctions) : numFuncs(numFunctions) {
    memoizedVal = 0;
    numPasses = numFuncs/5 + 1;
    memoizedHashes = gm_calloc<uint32_t>(numPasses*5);  // always > than multiple of buffers
}

uint64_t SHA1HashFamily::hash(uint32_t id, uint64_t val) {
    assert(id >= 0 && id < (uint32_t)numFuncs);
    if (val == memoizedVal) {
        //info("Memo hit 0x%x", memoizedHashes[id]);
        return (uint64_t) memoizedHashes[id];
    } else {
        uint64_t buffer[16];
        //sha1_context ctx;
        for (int i = 0; i < 16; i++) {
            buffer[i] = val;
        }

        for (int i = 0; i < numPasses; i++) {
            if (i > 0) { //change source
                for (int j = 0; j < 5; j++) {
                    buffer[j] ^= memoizedHashes[(i-1)*5 + j];
                }
            }
            sha1((unsigned char*) buffer, sizeof(buffer), (unsigned char*) &(memoizedHashes[i*5]));
        }
        /*info("SHA1: 0x%lx:", val);
        for (int i = 0; i < numFuncs; i++) {
            info(" %d: 0x%x", i, memoizedHashes[i]);
        }*/

        memoizedVal = val;
        return (uint64_t) memoizedHashes[id];
    }
}

#else  // _WITH_POLARSSL_

SHA1HashFamily::SHA1HashFamily(int numFunctions) {
    panic("Cannot use SHA1HashFamily, zsim was not compiled with PolarSSL");
}

uint64_t SHA1HashFamily::hash(uint32_t id, uint64_t val) {
    panic("???");
    return 0;
}

Qarma64HashFamily::Qarma64HashFamily(unsigned mask) {
    setMask = mask;
}

void Qarma64HashFamily::text2cell(cell_t* cell, uint64_t is) {
    // for 64 bits
    char* byte_ptr = (char*)&is;
    for (int i = 0; i < MAX_LENGTH / 8; i++) {
        char byte = byte_ptr[i];
        cell[2 * (7 - i) + 0] = (byte & 0xF0) >> 4;
        cell[2 * (7 - i) + 1] = byte & 0xF;
    }
}

uint64_t Qarma64HashFamily::cell2text(cell_t* cell) {
    uint64_t is = 0;
    for (int i = 0; i < MAX_LENGTH / 8; i++) {
        uint64_t byte = 0;
        byte = (cell[2 * i] << 4) | cell[2 * i + 1];
        is = is | (byte << (7 - i) * 8UL);
    }
    return is;
}

uint64_t Qarma64HashFamily::pseudo_reflect(uint64_t is, uint64_t tk) {
    cell_t cell[16];
    text2cell(cell, is);

    // ShuffleCells
    cell_t perm[16];
    for (int i = 0; i < 16; i++)
        perm[i] = cell[t[i]];

    // MixColumns
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            cell_t temp = 0;
            for (int j = 0; j < 4; j++) {
                int b;
                if ((b = M[4 * x + j]) != 0) {
                    cell_t a = perm[4 * j + y];
                    temp ^= ((a << b) & 0x0F) |
                        (a >> (4 - b));
                }
            }
            cell[4 * x + y] = temp;
        }
    }

    // AddRoundTweakey
    for (int i = 0; i < 16; i++)
        cell[i] ^= (tk >> (4 * (15 - i))) & 0xF;

    // ShuffleCells invert
    for (int i = 0; i < 16; i++)
        perm[i] = cell[t_inv[i]];

    return cell2text(perm);
}

uint64_t Qarma64HashFamily::forward(uint64_t is, uint64_t tk, int r) {
    is ^= tk;
    cell_t cell[16];
    text2cell(cell, is);

    if (r != 0) {
        // ShuffleCells
        cell_t perm[16];
        for (int i = 0; i < 16; i++)
            perm[i] = cell[t[i]];

        // MixColumns
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                cell_t temp = 0;
                for (int j = 0; j < 4; j++) {
                    int b;
                    if ((b = M[4 * x + j]) != 0) {
                        cell_t a = perm[4 * j + y];
                        temp ^= ((a << b) & 0x0F) |
                            (a >> (4 - b));
                    }
                }
                cell[4 * x + y] = temp;
            }
        }
    }

    // SubCells
    for (int i = 0; i < 16; i++) {
        cell[i] = subcells[cell[i]];
    }
    is = cell2text(cell);

    return is;
}

uint64_t Qarma64HashFamily::backward(uint64_t is, uint64_t tk, int r) {
    cell_t cell[16];
    text2cell(cell, is);

    // SubCells
    for (int i = 0; i < 16; i++) {
        cell[i] = subcells_inv[cell[i]];
    }

    if (r != 0) {
        cell_t mixc[16];
        // MixColumns
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                cell_t temp = 0;
                for (int j = 0; j < 4; j++) {
                    int b;
                    if ((b = M[4 * x + j]) != 0) {
                        cell_t a = cell[4 * j + y];
                        temp ^= ((a << b) & 0x0F) |
                            (a >> (4 - b));
                    }
                }
                mixc[4 * x + y] = temp;
            }
        }

        // ShuffleCells
        for (int i = 0; i < 16; i++)
            cell[i] = mixc[t_inv[i]];
    }

    is = cell2text(cell);
    is ^= tk;

    return is;
}

cell_t Qarma64HashFamily::LFSR(cell_t x) {
    cell_t b0 = (x >> 0) & 1;
    cell_t b1 = (x >> 1) & 1;
    cell_t b2 = (x >> 2) & 1;
    cell_t b3 = (x >> 3) & 1;

    return ((b0 ^ b1) << 3) | (b3 << 2) | (b2 << 1) | (b1 << 0);
}

cell_t Qarma64HashFamily::LFSR_inv(cell_t x) {
    cell_t b0 = (x >> 0) & 1;
    cell_t b1 = (x >> 1) & 1;
    cell_t b2 = (x >> 2) & 1;
    cell_t b3 = (x >> 3) & 1;

    return ((b0 ^ b3) << 0) | (b0 << 1) | (b1 << 2) | (b2 << 3);
}

uint64_t Qarma64HashFamily::forward_update_key(uint64_t T) {
    cell_t cell[16], temp[16];
    text2cell(cell, T);

    // h box
    for (int i = 0; i < 16; i++) {
        temp[i] = cell[h[i]];
    }

    // w LFSR
    temp[0] = LFSR(temp[0]);
    temp[1] = LFSR(temp[1]);
    temp[3] = LFSR(temp[3]);
    temp[4] = LFSR(temp[4]);
    temp[8] = LFSR(temp[8]);
    temp[11] = LFSR(temp[11]);
    temp[13] = LFSR(temp[13]);

    return cell2text(temp);
}

uint64_t Qarma64HashFamily::backward_update_key(uint64_t T) {
    cell_t cell[16], temp[16];
    text2cell(cell, T);

    // w LFSR invert
    cell[0] = LFSR_inv(cell[0]);
    cell[1] = LFSR_inv(cell[1]);
    cell[3] = LFSR_inv(cell[3]);
    cell[4] = LFSR_inv(cell[4]);
    cell[8] = LFSR_inv(cell[8]);
    cell[11] = LFSR_inv(cell[11]);
    cell[13] = LFSR_inv(cell[13]);

    // h box
    for (int i = 0; i < 16; i++) {
        temp[i] = cell[h_inv[i]];
    }

    return cell2text(temp);
}

uint64_t Qarma64HashFamily::qarma64_enc(uint64_t plaintext, uint64_t tweak, int rounds) {
    uint64_t w1 = ((w0 >> 1) | (w0 << (64 - 1))) ^ (w0 >> (16 * m - 1));
    uint64_t k1 = k0;

    uint64_t is = plaintext ^ w0;

    for (int i = 0; i < rounds; i++) {
        is = forward(is, k0 ^ tweak ^ c[i], i);
        tweak = forward_update_key(tweak);
    }

    is = forward(is, w1 ^ tweak, 1);
    is = pseudo_reflect(is, k1);
    is = backward(is, w0 ^ tweak, 1);

    for (int i = rounds - 1; i >= 0; i--) {
        tweak = backward_update_key(tweak);
        is = backward(is, k0 ^ tweak ^ c[i] ^ alpha, i);
    }

    is ^= w1;

    return is;
}

uint64_t Qarma64HashFamily::qarma64_dec(uint64_t plaintext, uint64_t tweak, int rounds) {
    uint64_t w1 = w0;
    w0 = ((w0 >> 1) | (w0 << (64 - 1))) ^ (w0 >> (16 * m - 1));

    cell_t k0_cell[16], k1_cell[16];
    text2cell(k0_cell, k0);
    // MixColumns
    for (int x = 0; x < 4; x++) {
        for (int y = 0; y < 4; y++) {
            cell_t temp = 0;
            for (int j = 0; j < 4; j++) {
                int b;
                if ((b = M[4 * x + j]) != 0) {
                    cell_t a = k0_cell[4 * j + y];
                    temp ^= ((a << b) & 0x0F) |
                        (a >> (4 - b));
                }
            }
            k1_cell[4 * x + y] = temp;
        }
    }
    uint64_t k1 = cell2text(k1_cell);

    k0 ^= alpha;

    uint64_t is = plaintext ^ w0;

    for (int i = 0; i < rounds; i++) {
        is = forward(is, k0 ^ tweak ^ c[i], i);
        tweak = forward_update_key(tweak);
    }

    is = forward(is, w1 ^ tweak, 1);
    is = pseudo_reflect(is, k1);
    is = backward(is, w0 ^ tweak, 1);

    for (int i = rounds - 1; i >= 0; i--) {
        tweak = backward_update_key(tweak);
        is = backward(is, k0 ^ tweak ^ c[i] ^ alpha, i);
    }

    is ^= w1;

    return is;
}

uint64_t Qarma64HashFamily::hash(uint32_t id, uint64_t val){
    // id: wayIdx, val: lineAddr
    /**
     * TODO: concatenate SDID (maybe also SecID) with way index as the tweak.
     */

    /**
     * SCv1: simply use both tag and index bits as the plaintext, whose output
     * has potential birthday-bound complexity.
     */
    // return qarma64_enc(val, id, NUM_ENC_ROUNDS);

    /**
     * SCv2: only use index bits as the plaintext, while the tag bits
     * constitute the tweak in order to mitigate birthday-bound index
     * collisions.
     */
    uint64_t indexBits = (val & setMask);
    uint64_t tweak = ((val & (~setMask)) | id);
    return qarma64_enc(indexBits, tweak, NUM_ENC_ROUNDS);
}

#endif  // _WITH_POLARSSL_
