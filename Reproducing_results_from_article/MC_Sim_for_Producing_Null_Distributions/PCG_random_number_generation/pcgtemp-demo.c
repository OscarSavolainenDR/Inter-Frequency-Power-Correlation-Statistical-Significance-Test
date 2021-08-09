/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

/*
 * This file was mechanically generated from tests/check-pcg32.c
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#include "pcg_basic.h"

int main(int argc, char** argv)
{
    // Read command-line options
    
    // 1st input is the length of the sequence
    // 2nd will be seed1
    // 3rd will be seed2
    // if not 2nd and 3rd inputs, sequence is deterministic

    int seq_length = 5;
    bool nondeterministic_seed = false;
    int i, seed_1, seed_2;
    
    if (argc>0) {
        printf("%d",argc);
        printf("\n");
    }

    //++argv;
    //--argc;
    //if (argc > 0 && strcmp(argv[0], "-r") == 0) {
    //    nondeterministic_seed = true;
        //++argv;
        //--argc;
    //}
    if (argc > 1) {
        seq_length = atoi(argv[1]);
        printf("%d",seq_length);
        printf("\n");
    }
    if (argc > 2) {
        nondeterministic_seed = true;
        seed_1 = atoi(argv[2]);
        seed_2 = atoi(argv[3]);
    }

    // In this version of the code, we'll use a local rng, rather than the
    // global one.

    pcg32_random_t rng;

    // You should *always* seed the RNG.  The usual time to do it is the
    // point in time when you create RNG (typically at the beginning of the
    // program).
    //
    // pcg32_srandom_r takes two 64-bit constants (the initial state, and the
    // rng sequence selector; rngs with different sequence selectors will
    // *never* have random sequences that coincide, at all) - the code below
    // shows three possible ways to do so.

    if (nondeterministic_seed) {
        // Seed with external entropy -- the time and some program addresses
        // (which will actually be somewhat random on most modern systems).
        // A better solution, entropy_getbytes, using /dev/random, is provided
        // in the full library.
        
        //pcg32_srandom_r(&rng, seed_1, (intptr_t)&seq_length);
        pcg32_srandom_r(&rng, seed_1, seed_2);
                                         //pcg32_srandom_r(&rng, seed_1, seed_2);
    } else {
        // Seed with a fixed constant

        // pcg32_srandom_r(&rng, 42u, 54u);
        pcg32_srandom_r(&rng, 42u, 54u);
    }


    printf("  32bit:");
    printf("\n");
    for (i = 1; i <= seq_length; ++i) {
        //nb = (int)strtol(pcg32_random_r(&rng), NULL, 0);
        printf(" %d", pcg32_random_r(&rng));
        //printf("%*c", 2, ' ');
        printf(" | ");
       // printf("\n");
        
    }
    printf("Done \n");

    return 0;
}
