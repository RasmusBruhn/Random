#include <stdio.h>
#include <stdint.h>
#include "Random.h"
#include <time.h>
#include <math.h>

typedef struct __NormalParams 
{
    double mu;
    double sigma;
} NormalParams;

void PrintIntArray(char *Name, uint64_t *Array, size_t Length);

void PrintFloatArray(char *Name, double *Array, size_t Length);

double *HistInt(uint64_t Min, uint64_t Max, uint32_t Bins, uint64_t *Array, size_t Size);

double *HistFloat(double Min, double Max, uint32_t Bins, double *Array, size_t Size);

double *linspace(double Min, double Max, uint32_t Length);

uint64_t *arange(uint64_t Min, uint64_t Step, uint32_t Length);

double NormalPDF(double x, void *Params);

double ExpPDF(double x, void *Params);

double ExpSampling(RNG_Seed Seed, void *Params);

double NormalSampling(RNG_Seed Seed, double Mu, double Sigma);

void NormalPDF_ArrayM(double *x, void *Params, double *Array, size_t Size);

void ExpPDF_ArrayM(double *x, void *Params, double *Array, size_t Size);

void ExpSampling_ArrayM(RNG_Seed Seed, void *Params, double *Array, size_t Size);

double *NormalSampling_Array(RNG_Seed Seed, double Mu, double Sigma, size_t Size);

// Test Random.h
int main(int argc, char **argv)
{
    size_t Size = 1000000;
    double *Hist = NULL;
    uint64_t *IntArray = NULL;
    double *FloatArray = NULL;
    uint32_t Bins = 11;
    double *xFloat = NULL;
    uint64_t *xInt = NULL;

    // Generate seed
    RNG_Seed Seed = RNG_GenSeed();
    printf("Seed: %lu\n", *Seed);

    // Generate uniform ints
    printf("Uniform ints: %lu, %lu, %lu\n", RNG_Int(Seed, 0, Bins - 1), RNG_Int(Seed, 0, Bins - 1), RNG_Int(NULL, 0, Bins - 1));

    IntArray = RNG_Int_Array(Seed, 0, Bins - 1, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Uniform int distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    // Generate uniform floats
    printf("Uniform floats: %.3g, %.3g, %.3g\n", RNG_Float(Seed, 0, 1), RNG_Float(Seed, 0, 1), RNG_Float(NULL, 0, 1));

    FloatArray = RNG_Float_Array(Seed, 0, 1, Size);
    Hist = HistFloat(0, 1, Bins, FloatArray, Size);
    PrintFloatArray("Uniform float distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    // Generate exp
    printf("Exp: %.3g, %.3g, %.3g\n", RNG_Exp(Seed, 0, 1), RNG_Exp(Seed, 0, 1), RNG_Exp(NULL, 0, 1));

    FloatArray = RNG_Exp_Array(Seed, 0, 1, Size);
    Hist = HistFloat(0, 5, Bins, FloatArray, Size);
    PrintFloatArray("Exp distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    printf("Exp PDF: %.3g, %.3g, %.3g\n", RNG_ExpPDF(0, 0, 1), RNG_ExpPDF(1, 0, 1), RNG_ExpPDF(2, 0, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpPDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp PDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp CDF: %.3g, %.3g, %.3g\n", RNG_ExpCDF(0, 0, 1), RNG_ExpCDF(1, 0, 1), RNG_ExpCDF(2, 0, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpCDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp CDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp ICDF: %.3g, %.3g, %.3g\n", RNG_ExpICDF(0, 0, 1), RNG_ExpICDF(0.4, 0, 1), RNG_ExpICDF(0.8, 0, 1));

    xFloat = linspace(0, 1, Bins);
    FloatArray = RNG_ExpICDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp ICDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    // Generate normal
    printf("Normal: %.3g, %.3g, %.3g\n", RNG_Normal(Seed, 0, 1), RNG_Normal(Seed, 0, 1), RNG_Normal(NULL, 0, 1));

    FloatArray = RNG_Normal_Array(Seed, 0, 1, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Normal distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    printf("Normal PDF: %.3g, %.3g, %.3g\n", RNG_NormalPDF(-1, 0, 1), RNG_NormalPDF(0, 0, 1), RNG_NormalPDF(1, 0, 1));

    xFloat = linspace(-3, 3, Bins);
    FloatArray = RNG_NormalPDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal PDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Normal CDF: %.3g, %.3g, %.3g\n", RNG_NormalCDF(-1, 0, 1), RNG_NormalCDF(0, 0, 1), RNG_NormalCDF(1, 0, 1));

    xFloat = linspace(-3, 3, Bins);
    FloatArray = RNG_NormalCDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal CDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Normal ICDF: %.3g, %.3g, %.3g\n", RNG_NormalICDF(0.1, 0, 1), RNG_NormalICDF(0.5, 0, 1), RNG_NormalICDF(0.9, 0, 1));

    xFloat = linspace(0, 1, Bins);
    FloatArray = RNG_NormalICDF_Array(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal ICDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    // Generate normal using monte carlo
    printf("Normal: %.3g, %.3g, %.3g\n", NormalSampling(Seed, 0, 1), NormalSampling(Seed, 0, 1), NormalSampling(NULL, 0, 1));

    FloatArray = NormalSampling_Array(Seed, 0, 1, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Normal distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    NormalParams Params = {.mu = 0, .sigma = 1};
    FloatArray = malloc(sizeof(double) * Size);
    ExpSampling_ArrayM(Seed, &Params, FloatArray, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Double exp distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    xFloat = linspace(-3, 3, Bins);
    ExpPDF_ArrayM(xFloat, &Params, xFloat, Bins);
    PrintFloatArray("Double exp PDF distribution", xFloat, Bins);

    free(xFloat);
    free(FloatArray);

    RNG_DestroySeed(Seed);

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

double NormalPDF(double x, void *Params)
{
    return RNG_NormalPDF(x, ((NormalParams *)Params)->mu, ((NormalParams *)Params)->sigma);
}

double ExpPDF(double x, void *Params)
{
    return M_2_SQRTPI / (((NormalParams *)Params)->sigma * 2 * M_SQRT2) * exp(-fabs((x - ((NormalParams *)Params)->mu) / ((NormalParams *)Params)->sigma) + 0.5);
}

double ExpSampling(RNG_Seed Seed, void *Params)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FastFloat(Seed);

    double Sign = 1;

    if (Uniform >= 0.5)
    {
        Sign = -1;
        Uniform -= 0.5;
    }

    // Get from current distribution
    return ((NormalParams *)Params)->mu - Sign * ((NormalParams *)Params)->sigma * log(1 - 2 * Uniform);
}

double NormalSampling(RNG_Seed Seed, double Mu, double Sigma)
{
    NormalParams Params = {.mu = Mu, .sigma = Sigma};

    return RNG_MonteCarlo(Seed, &NormalPDF, &ExpPDF, &ExpSampling, &Params);
}

void NormalPDF_ArrayM(double *x, void *Params, double *Array, size_t Size)
{
    RNG_NormalPDF_ArrayM(x, ((NormalParams *)Params)->mu, ((NormalParams *)Params)->sigma, Array, Size);
}

void ExpPDF_ArrayM(double *x, void *Params, double *Array, size_t Size)
{
    double A = 1 / (((NormalParams *)Params)->sigma * sqrt(2 * M_PI)) * exp(0.5);
    double B = 1 / ((NormalParams *)Params)->sigma;

    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-fabs(B * (*x - ((NormalParams *)Params)->mu)));
}

void ExpSampling_ArrayM(RNG_Seed Seed, void *Params, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FastFloat(Seed);

        double Sign = ((NormalParams *)Params)->sigma;

        if (Uniform >= 0.5)
        {
            Sign *= -1;
            Uniform -= 0.5;
        }

        // Get from distribution
        *List = ((NormalParams *)Params)->mu - Sign * log(1 - 2 * Uniform);
    }
}

double *NormalSampling_Array(RNG_Seed Seed, double Mu, double Sigma, size_t Size)
{
    NormalParams Params = {.mu = Mu, .sigma = Sigma};

    return RNG_MonteCarlo_Array(Seed, &NormalPDF_ArrayM, &ExpPDF_ArrayM, &ExpSampling_ArrayM, &Params, Size);
}
