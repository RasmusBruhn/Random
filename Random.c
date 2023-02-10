#include <stdio.h>
#include <stdint.h>
#include "Random.h"
#include <time.h>

// Test Random.h
int main(int argc, char **argv)
{
    // Generate seed
    uint64_t Seed = RNG_GenSeed();
    printf("Seed: %u\n", Seed);

    // Generate some random integers
    printf("Random ints: %X, %X, %X\n", RNG_Int(&Seed), RNG_Int(&Seed), RNG_Int(NULL));

    // Generate some random floats
    printf("Random floats: %.3g, %.3g, %.3g\n", RNG_Float(&Seed), RNG_Float(&Seed), RNG_Float(NULL));

    size_t Size = 10;

    // Generate int list
    uint64_t *IntArray = RNG_Int_Array(&Seed, Size);

    printf("Random int array: [");

    for (uint64_t *List = IntArray, *EndList = IntArray + Size - 1; List < EndList; ++List)
        printf("%u, ", *List % 100);

    printf("%u]\n", IntArray[Size - 1] % 100);

    // Generate int list
    double *FloatArray = RNG_Float_Array(&Seed, Size);

    printf("Random float array: [");

    for (double *List = FloatArray, *EndList = FloatArray + Size - 1; List < EndList; ++List)
        printf("%.3g, ", *List);

    printf("%.3g]\n", FloatArray[Size - 1]);

    free(IntArray);
    free(FloatArray);

    return 0;
}