#ifndef RANDOM2_H_INCLUDED
#define RANDOM2_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

// Type to store seed in
typedef uint64_t RNG_Seed;

// Constants
#define RNG_MOD 0x7FFFFFFF
#define RNG_SIZE 0x80000000
#define RNG_MULTIPLIER 0x7FFFFFED
#define RNG_CONSTANT 0x7FFFFFC3

uint64_t RNG_GlobalSeed = 0;

// Get a random uint32_t and update the seed to be that number
#define RNG_RandS(Seed) ((uint32_t) ((Seed) = ((Seed) * RNG_MULTIPLIER + RNG_CONSTANT) % RNG_MOD))

// Get a random  double between 0 inclusive and 1 exclusive and update seed
#define RNG_RandSf(Seed) ((double) ((Seed) = ((Seed) * RNG_MULTIPLIER + RNG_CONSTANT) % RNG_MOD) / (double) RNG_SIZE)

// Initialize global random numbers
#define RNG_Init() extern RNG_Seed RNG_GlobalSeed

// Set global seed
#define RNG_SetSeed(Seed) RNG_GlobalSeed = (Seed)

// Set global seed to random seed
#define RNG_RandSeed() RNG_SetSeed(time(NULL))

// Get the global seed
#define RNG_GetSeed() RNG_GlobalSeed

// Get a random uint32_t
#define RNG_Rand() RNG_RandS(RNG_GlobalSeed)

// Get a random double between 0 and 1 exclusive
#define RNG_Randf() RNG_RandSf(RNG_GlobalSeed)

#endif // RANDOM2_H_INCLUDED
