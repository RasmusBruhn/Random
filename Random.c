#include <stdio.h>
#include <stdint.h>
#include "Random.h"
#include <time.h>

void PrintIntArray(char *Name, uint64_t *Array, size_t Length);

void PrintFloatArray(char *Name, double *Array, size_t Length);

double *HistInt(uint64_t Min, uint64_t Max, uint32_t Bins, uint64_t *Array, size_t Size);

double *HistFloat(double Min, double Max, uint32_t Bins, double *Array, size_t Size);

double *linspace(double Min, double Max, uint32_t Length);

uint64_t *arange(uint64_t Min, uint64_t Step, uint32_t Length);

// Test Random.h
int main(int argc, char **argv)
{
    size_t Size = 100000;
    double *Hist = NULL;
    uint64_t *IntArray = NULL;
    double *FloatArray = NULL;
    uint32_t Bins = 10;
    double *xFloat = NULL;
    uint64_t *xInt = NULL;

    // Generate seed
    uint64_t Seed = RNG_GenSeed();
    printf("Seed: %lu\n", Seed);

    // Generate uniform ints
    printf("Uniform ints: %lu, %lu, %lu\n", RNG_Int(&Seed, 0, Bins - 1), RNG_Int(&Seed, 0, Bins - 1), RNG_Int(NULL, 0, Bins - 1));

    IntArray = RNG_Int_Array(&Seed, 0, Bins - 1, Size);
    Hist = HistInt(0, 9, Bins, IntArray, Size);
    PrintFloatArray("Uniform int distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    // Generate uniform floats
    printf("Uniform floats: %.3g, %.3g, %.3g\n", RNG_Float(&Seed, 0, 1), RNG_Float(&Seed, 0, 1), RNG_Float(NULL, 0, 1));

    FloatArray = RNG_Float_Array(&Seed, 0, 1, Size);
    Hist = HistFloat(0, 1, Bins, FloatArray, Size);
    PrintFloatArray("Uniform float distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    // Generate exp
    printf("Exp: %.3g, %.3g, %.3g\n", RNG_Exp(&Seed, 1), RNG_Exp(&Seed, 1), RNG_Exp(NULL, 1));

    FloatArray = RNG_Exp_Array(&Seed, 1, Size);
    Hist = HistFloat(0, 5, Bins, FloatArray, Size);
    PrintFloatArray("Exp distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    printf("Exp PDF: %.3g, %.3g, %.3g\n", RNG_ExpPDF(0, 1), RNG_ExpPDF(1, 1), RNG_ExpPDF(2, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpPDF_Array(xFloat, 1, Size);
    PrintFloatArray("Exp PDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp CDF: %.3g, %.3g, %.3g\n", RNG_ExpCDF(0, 1), RNG_ExpCDF(1, 1), RNG_ExpCDF(2, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpCDF_Array(xFloat, 1, Size);
    PrintFloatArray("Exp CDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp ICDF: %.3g, %.3g, %.3g\n", RNG_ExpICDF(0, 1), RNG_ExpICDF(0.4, 1), RNG_ExpICDF(0.8, 1));

    xFloat = linspace(0, 1, Bins);
    FloatArray = RNG_ExpICDF_Array(xFloat, 1, Size);
    PrintFloatArray("Exp ICDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    return 0;
}

void PrintIntArray(char *Name, uint64_t *Array, size_t Length)
{
    printf("%s: [", Name);

    for (uint64_t *List = Array, *EndList = Array + Length - 1; List < EndList; ++List)
        printf("%lu, ", *List);

    printf("%lu]\n", Array[Length - 1]);
}

void PrintFloatArray(char *Name, double *Array, size_t Length)
{
    printf("%s: [", Name);

    for (double *List = Array, *EndList = Array + Length - 1; List < EndList; ++List)
        printf("%.3g, ", *List);

    printf("%.3g]\n", Array[Length - 1]);
}

double *HistInt(uint64_t Min, uint64_t Max, uint32_t Bins, uint64_t *Array, size_t Size)
{
    // Get bins
    double *BinCounts = (double *)malloc(sizeof(double) * Bins);

    for (double *List = BinCounts, *EndList = BinCounts + Bins; List < EndList; ++List)
        *List = 0;

    // Fill them up
    for (uint64_t *Numbers = Array, *EndNumbers = Array + Size; Numbers < EndNumbers; ++Numbers) 
    {
        int64_t Index = (*Numbers - Min) * Bins / (1 + Max - Min);

        if (Index >= 0 && Index < Bins)
            BinCounts[Index] += 1;
    }

    for (double *List = BinCounts, *EndList = BinCounts + Bins; List < EndList; ++List)
        *List /= Size * (1 + Max - Min) / Bins;

    return BinCounts;
}

double *HistFloat(double Min, double Max, uint32_t Bins, double *Array, size_t Size)
{
    // Get bins
    double *BinCounts = (double *)malloc(sizeof(double) * Bins);

    for (double *List = BinCounts, *EndList = BinCounts + Bins; List < EndList; ++List)
        *List = 0;

    // Fill them up
    for (double *Numbers = Array, *EndNumbers = Array + Size; Numbers < EndNumbers; ++Numbers)
    {
        int64_t Index = (uint64_t) ((*Numbers - Min) * Bins / (Max - Min));

        if (Index >= 0 && Index < Bins)
            BinCounts[Index] += 1;
    }

    for (double *List = BinCounts, *EndList = BinCounts + Bins; List < EndList; ++List)
        *List /= Size * (Max - Min) / Bins;

    return BinCounts;
}

double *linspace(double Min, double Max, uint32_t Length)
{
    double *Values = (double *)malloc(sizeof(double) * Length);

    double Step = (Max - Min) / (double) Length;
    Min += Step / 2;
    Max -= Step / 2;
    double Prev = Min - Step;

    for (double *List = Values, *EndList = Values + Length; List < EndList; ++List)
    {
        Prev += Step;
        *List = Prev;
    }

    return Values;
}

uint64_t *arange(uint64_t Min, uint64_t Step, uint32_t Length)
{
    uint64_t *Values = (uint64_t *)malloc(sizeof(uint64_t) * Length);

    uint64_t Prev = Min - Step;

    for (uint64_t *List = Values, *EndList = Values + Length; List < EndList; ++List)
    {
        Prev += Step;
        *List = Prev;
    }

    return Values;
}
