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

void NormalPDFArrayM(double *x, void *Params, double *Array, size_t Size);

void ExpPDFArrayM(double *x, void *Params, double *Array, size_t Size);

void ExpSamplingArrayM(RNG_Seed Seed, void *Params, double *Array, size_t Size);

double *NormalSamplingArray(RNG_Seed Seed, double Mu, double Sigma, size_t Size);

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
    RNG_Seed Seed = RNG_SeedGenerate();
    printf("Seed: %lu\n", *Seed);

    // Generate uniform ints
    printf("Uniform ints: %lu, %lu, %lu\n", RNG_Int(Seed, 0, Bins - 1), RNG_Int(Seed, 0, Bins - 1), RNG_Int(NULL, 0, Bins - 1));

    IntArray = RNG_IntArray(Seed, 0, Bins - 1, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Uniform int distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    // Generate uniform floats
    printf("Uniform floats: %.3g, %.3g, %.3g\n", RNG_Float(Seed, 0, 1), RNG_Float(Seed, 0, 1), RNG_Float(NULL, 0, 1));

    FloatArray = RNG_FloatArray(Seed, 0, 1, Size);
    Hist = HistFloat(0, 1, Bins, FloatArray, Size);
    PrintFloatArray("Uniform float distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    // Generate exp
    printf("Exp: %.3g, %.3g, %.3g\n", RNG_Exp(Seed, 0, 1), RNG_Exp(Seed, 0, 1), RNG_Exp(NULL, 0, 1));

    FloatArray = RNG_ExpArray(Seed, 0, 1, Size);
    Hist = HistFloat(0, 5, Bins, FloatArray, Size);
    PrintFloatArray("Exp distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    printf("Exp PDF: %.3g, %.3g, %.3g\n", RNG_ExpPDF(0, 0, 1), RNG_ExpPDF(1, 0, 1), RNG_ExpPDF(2, 0, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpPDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp PDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp CDF: %.3g, %.3g, %.3g\n", RNG_ExpCDF(0, 0, 1), RNG_ExpCDF(1, 0, 1), RNG_ExpCDF(2, 0, 1));

    xFloat = linspace(0, 5, Bins);
    FloatArray = RNG_ExpCDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp CDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Exp ICDF: %.3g, %.3g, %.3g\n", RNG_ExpICDF(0, 0, 1), RNG_ExpICDF(0.4, 0, 1), RNG_ExpICDF(0.8, 0, 1));

    xFloat = linspace(0, 1, Bins);
    FloatArray = RNG_ExpICDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Exp ICDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    // Generate normal
    printf("Normal: %.3g, %.3g, %.3g\n", RNG_Normal(Seed, 0, 1), RNG_Normal(Seed, 0, 1), RNG_Normal(NULL, 0, 1));

    FloatArray = RNG_NormalArray(Seed, 0, 1, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Normal distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    printf("Normal PDF: %.3g, %.3g, %.3g\n", RNG_NormalPDF(-1, 0, 1), RNG_NormalPDF(0, 0, 1), RNG_NormalPDF(1, 0, 1));

    xFloat = linspace(-3, 3, Bins);
    FloatArray = RNG_NormalPDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal PDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Normal CDF: %.3g, %.3g, %.3g\n", RNG_NormalCDF(-1, 0, 1), RNG_NormalCDF(0, 0, 1), RNG_NormalCDF(1, 0, 1));

    xFloat = linspace(-3, 3, Bins);
    FloatArray = RNG_NormalCDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal CDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    printf("Normal ICDF: %.3g, %.3g, %.3g\n", RNG_NormalICDF(0.1, 0, 1), RNG_NormalICDF(0.5, 0, 1), RNG_NormalICDF(0.9, 0, 1));

    xFloat = linspace(0, 1, Bins);
    FloatArray = RNG_NormalICDFArray(xFloat, 0, 1, Bins);
    PrintFloatArray("Normal ICDF distribution", FloatArray, Bins);

    free(xFloat);
    free(FloatArray);

    // Generate normal using monte carlo
    printf("Normal: %.3g, %.3g, %.3g\n", NormalSampling(Seed, 0, 1), NormalSampling(Seed, 0, 1), NormalSampling(NULL, 0, 1));

    FloatArray = NormalSamplingArray(Seed, 0, 1, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Normal distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    NormalParams Params = {.mu = 0, .sigma = 1};
    FloatArray = malloc(sizeof(double) * Size);
    ExpSamplingArrayM(Seed, &Params, FloatArray, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Double exp distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    xFloat = linspace(-3, 3, Bins);
    ExpPDFArrayM(xFloat, &Params, xFloat, Bins);
    PrintFloatArray("Double exp PDF distribution", xFloat, Bins);

    free(xFloat);

    RNG_SeedDestroy(Seed);

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
    return 1 / (((NormalParams *)Params)->sigma * _RNG_SQRTPI * _RNG_SQRT2) * exp(-fabs((x - ((NormalParams *)Params)->mu) / ((NormalParams *)Params)->sigma) + 0.5);
}

double ExpSampling(RNG_Seed Seed, void *Params)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FloatFast(Seed);

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

void NormalPDFArrayM(double *x, void *Params, double *Array, size_t Size)
{
    RNG_NormalPDFArrayM(x, ((NormalParams *)Params)->mu, ((NormalParams *)Params)->sigma, Array, Size);
}

void ExpPDFArrayM(double *x, void *Params, double *Array, size_t Size)
{
    double A = 1 / (((NormalParams *)Params)->sigma * sqrt(2 * 3.141592653589793)) * exp(0.5);
    double B = 1 / ((NormalParams *)Params)->sigma;

    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-fabs(B * (*x - ((NormalParams *)Params)->mu)));
}

void ExpSamplingArrayM(RNG_Seed Seed, void *Params, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FloatFast(Seed);

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

double *NormalSamplingArray(RNG_Seed Seed, double Mu, double Sigma, size_t Size)
{
    NormalParams Params = {.mu = Mu, .sigma = Sigma};

    return RNG_MonteCarloArray(Seed, &NormalPDFArrayM, &ExpPDFArrayM, &ExpSamplingArrayM, &Params, Size);
}
