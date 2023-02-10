#ifndef RANDOM2_H_INCLUDED
#define RANDOM2_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define ERR_PREFIX RNG
#include <Error.h>

enum _RNG_ErrorID {
    _RNG_ERRORID_GENERAL_NONE = 0x700000000,
    _RNG_ERRORID_INTARRAY_MALLOC = 0x700010200,
    _RNG_ERRORID_FLOATARRAY_MALLOC = 0x700020200
};

#define _RNG_ERRORMES_MALLOC "Unable to allocate memory (Size: %lu)"

// Constants
#define RNG_MAX 18446744073709551616.
#define RNG_MULTIPLIER 6364136223846793005
#define RNG_CONSTANT 1442695040888963407

uint64_t _RNG_GlobalSeed = 0;

// Get a random uint32_t and update the seed to be that number
// Seed (uint64_t): The seed to use and update
#define RNG_FastInt(Seed) ((Seed) = (Seed) * RNG_MULTIPLIER + RNG_CONSTANT)

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (uint64_t): The seed to use and update
#define RNG_FastFloat(Seed) ((double) RNG_FastInt(Seed) / RNG_MAX)

// Generates a seed, returns it as a uint32
#define RNG_GenSeed() ((uint64_t) time(NULL) ^ (uint64_t) clock())


// Get a random uint64_t and update the seed to be that number
// Seed (uint64_t): The seed to use and update, NULL to use global seed
uint64_t RNG_Int(uint64_t *Seed)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return RNG_FastInt(*Seed);
}

// Get an array of random uint64_t and updates the seed
// Seed (uint64_t): The seed to use and update, NULL to use global seed
// Size (size_t): The size of the array
uint64_t *RNG_Int_Array(uint64_t *Seed, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_INTARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    // Fill memory
    for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = RNG_FastInt(*Seed);

    return Array;
}

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (uint64_t): The seed to use and update, NULL to use global seed
double RNG_Float(uint64_t *Seed)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return RNG_FastFloat(*Seed);
}

// Get an array of random double between 0 inclusive and 1 exclusive and update seed
// Seed (uint64_t): The seed to use and update, NULL to use global seed
// Size (size_t): The size of the array
double *RNG_Float_Array(uint64_t *Seed, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_FLOATARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = RNG_FastFloat(*Seed);

    return Array;
}

#endif // RANDOM2_H_INCLUDED
