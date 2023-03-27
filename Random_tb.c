#include <Debug2.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "Random.h"
#include <Defines.h>

#define _RNG_SQRT2 1.4142135623730951
#define _RNG_SQRTPI 1.7724538509055159
#define _RNG_1_E 0.36787944117144233
#define _RNG_E 2.718281828459045

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

double NormalPDF(double x, const void *Params);

double ExpPDF(double x, const void *Params);

double ExpSampling(RNG_Seed *Seed, const void *Params);

double NormalSampling(RNG_Seed *Seed, double Mu, double Sigma);

double *NormalPDFArrayM(double *x, const void *Params, double *Array, size_t Size);

double *ExpPDFArrayM(double *x, const void *Params, double *Array, size_t Size);

double *ExpSamplingArrayM(RNG_Seed *Seed, const void *Params, double *Array, size_t Size);

double *NormalSamplingArray(RNG_Seed *Seed, double Mu, double Sigma, size_t Size);

// Test Random.h
int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv)
{
    size_t Size = 1000000;
    size_t TestSize = 100;
    double *Hist = NULL;
    uint64_t *IntArray = NULL;
    double *FloatArray = NULL;
    uint32_t Bins = 11;
    double *xFloat = NULL;
    uint64_t *xInt = NULL;

    // Generate seed
    RNG_Seed *Seed = RNG_SeedGenerate();
    printf("Seed: " PRINT_UINT64 "\n", *Seed);

    // Generate uniform ints
    printf("Uniform ints: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Int(Seed, 0, Bins - 1), RNG_Int(Seed, 0, Bins - 1), RNG_Int(NULL, 0, Bins - 1));

    IntArray = RNG_IntArray(Seed, 0, Bins - 1, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Uniform int distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    // Test special cases for ints
    IntArray = RNG_IntArray(Seed, 0, 0, TestSize);
    for (uint64_t *List = IntArray, *EndList = IntArray + TestSize; List < EndList; ++List)
        if (*List != 0)
        {
            printf("Error generating span 0 uniform ints: " PRINT_UINT64 "\n", *List);
            return -1;
        }

    free(IntArray);

    // Test errors for ints
    if (RNG_Int(Seed, 10, 9) != (uint64_t)-1)
    {
        printf("Excepted error for negative span int");
        return -1;
    }

    if (RNG_IntArray(Seed, 10, 9, TestSize) != NULL)
    {
        printf("Excepted error for negative span int array");
        return -1;
    }

    RNG_ErrorArchive();
    RNG_ErrorArchive();

    // Generate uniform floats
    printf("Uniform floats: %.3g, %.3g, %.3g\n", RNG_Float(Seed, 0, 1), RNG_Float(Seed, 0, 1), RNG_Float(NULL, 0, 1));

    FloatArray = RNG_FloatArray(Seed, 0, 1, Size);
    Hist = HistFloat(0, 1, Bins, FloatArray, Size);
    PrintFloatArray("Uniform float distribution", Hist, Bins);

    free(FloatArray);
    free(Hist);

    // Test special cases for floats
    FloatArray = RNG_FloatArray(Seed, 0, 0, TestSize);
    for (double *List = FloatArray, *EndList = FloatArray + TestSize; List < EndList; ++List)
        if (*List != 0)
        {
            printf("Error generating span 0 uniform floats: %.3g\n", *List);
            return -1;
        }

    free(FloatArray);

    // Test errors for floats
    if (!isnan(RNG_Float(Seed, 1, 0)))
    {
        printf("Excepted error for negative span float");
        return -1;
    }

    if (RNG_FloatArray(Seed, 1, 0, TestSize) != NULL)
    {
        printf("Excepted error for negative span float array");
        return -1;
    }

    RNG_ErrorArchive();
    RNG_ErrorArchive();

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

    // Test special cases for exp
    FloatArray = RNG_ExpArray(Seed, 0, 0, TestSize);
    for (double *List = FloatArray, *EndList = FloatArray + TestSize; List < EndList; ++List)
        if (*List != 0)
        {
            printf("Error generating l = 0, exp: %.3g\n", *List);
            return -1;
        }

    free(FloatArray);

    if (RNG_ExpPDF(-1, 0, 1) != 0 || RNG_ExpPDF(INFINITY, 0, 1) != 0)
    {
        printf("Error evaluating PDF for exp");
        return -1;
    }

    if (RNG_ExpPDF(-1, 0, 0) != 0 || !isinf(RNG_ExpPDF(0, 0, 0)) || RNG_ExpPDF(1, 0, 0) != 0)
    {
        printf("Error evaluating PDF for l = 0, exp");
        return -1;
    }

    if (RNG_ExpCDF(-1, 0, 1) != 0 || RNG_ExpCDF(INFINITY, 0, 1) != 1)
    {
        printf("Error evaluating CDF for exp");
        return -1;
    }

    if (RNG_ExpCDF(-1, 0, 0) != 0 || RNG_ExpCDF(0, 0, 0) != 1 || RNG_ExpCDF(1, 0, 0) != 1)
    {
        printf("Error evaluating CDF for l = 0, exp");
        return -1;
    }

    if (RNG_ExpICDF(-1, 0, 1) != 0 || !isinf(RNG_ExpICDF(2, 0, 1)))
    {
        printf("Error evaluating ICDF for exp");
        return -1;
    }

    if (RNG_ExpICDF(-1, 0, 0) != 0 || RNG_ExpICDF(0, 0, 0) != 0 || RNG_ExpICDF(1, 0, 0) != 0)
    {
        printf("Error evaluating ICDF for l = 0, exp");
        return -1;
    }

    // Test errors for exp
    if (!isnan(RNG_Exp(Seed, 0, -1)))
    {
        printf("Excepted error for negative l exp");
        return -1;
    }

    if (RNG_ExpArray(Seed, 0, -1, TestSize) != NULL)
    {
        printf("Excepted error for negative l exp array");
        return -1;
    }

    if (!isnan(RNG_ExpPDF(0, 0, -1)))
    {
        printf("Excepted error for negative l exp PDF");
        return -1;
    }

    if (RNG_ExpPDFArray(0, 0, -1, TestSize) != NULL)
    {
        printf("Excepted error for negative l exp PDF array");
        return -1;
    }

    if (!isnan(RNG_ExpCDF(0, 0, -1)))
    {
        printf("Excepted error for negative l exp CDF");
        return -1;
    }

    if (RNG_ExpCDFArray(0, 0, -1, TestSize) != NULL)
    {
        printf("Excepted error for negative l exp CDF array");
        return -1;
    }

    if (!isnan(RNG_ExpICDF(0, 0, -1)))
    {
        printf("Excepted error for negative l exp ICDF");
        return -1;
    }

    if (RNG_ExpICDFArray(0, 0, -1, TestSize) != NULL)
    {
        printf("Excepted error for negative l exp ICDF array");
        return -1;
    }

    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();
    RNG_ErrorArchive();

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
    printf("Monte Carlo Normal: %.3g, %.3g, %.3g\n", NormalSampling(Seed, 0, 1), NormalSampling(Seed, 0, 1), NormalSampling(NULL, 0, 1));

    FloatArray = NormalSamplingArray(Seed, 0, 1, Size);
    Hist = HistFloat(-3, 3, Bins, FloatArray, Size);
    PrintFloatArray("Monte Carlo Normal distribution", Hist, Bins);

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

    // Generate Poisson
    printf("Poisson Large: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Poisson(Seed, 20), RNG_Poisson(Seed, 20), RNG_Poisson(NULL, 20));

    IntArray = RNG_PoissonArray(Seed, 20, Size);
    Hist = HistInt(20 - (Bins - 1) / 2, 20 + (Bins - 1) / 2, Bins, IntArray, Size);
    PrintFloatArray("Poisson Large distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    printf("Poisson Large PMF: %.3g, %.3g, %.3g\n", RNG_PoissonPMF(15, 20), RNG_PoissonPMF(20, 20), RNG_PoissonPMF(25, 20));

    xInt = arange(20 - (Bins - 1) / 2, 1, Bins);
    FloatArray = RNG_PoissonPMFArray(xInt, 20, Bins);
    PrintFloatArray("Poisson Large PMF distribution", FloatArray, Bins);

    free(xInt);
    free(FloatArray);

    printf("Poisson Small: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Poisson(Seed, 0.5), RNG_Poisson(Seed, 0.5), RNG_Poisson(NULL, 0.5));

    IntArray = RNG_PoissonArray(Seed, 0.5, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Poisson Small distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    printf("Poisson Small PMF: %.3g, %.3g, %.3g\n", RNG_PoissonPMF(0, 0.5), RNG_PoissonPMF(1, 0.5), RNG_PoissonPMF(2, 0.5));

    xInt = arange(0, 1, Bins);
    FloatArray = RNG_PoissonPMFArray(xInt, 0.5, Bins);
    PrintFloatArray("Poisson Small PMF distribution", FloatArray, Bins);

    free(xInt);
    free(FloatArray);

    printf("Binomial Normal: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Binomial(Seed, Bins - 1, 0.5), RNG_Binomial(Seed, Bins - 1, 0.5), RNG_Binomial(NULL, Bins - 1, 0.5));

    IntArray = RNG_BinomialArray(Seed, Bins - 1, 0.3, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Binomial Normal distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    printf("Binomial Normal PMF: %.3g, %.3g, %.3g\n", RNG_BinomialPMF(0, Bins - 1, 0.5), RNG_BinomialPMF((Bins - 1) / 2, Bins - 1, 0.5), RNG_BinomialPMF(Bins - 1, Bins - 1, 0.5));

    xInt = arange(0, 1, Bins);
    FloatArray = RNG_BinomialPMFArray(xInt, Bins - 1, 0.3, Bins);
    PrintFloatArray("Binomial Normal PMF distribution", FloatArray, Bins);

    free(xInt);
    free(FloatArray);

    printf("Binomial Small: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Binomial(Seed, Bins - 1, 0.02), RNG_Binomial(Seed, Bins - 1, 0.5), RNG_Binomial(NULL, Bins - 1, 0.02));

    IntArray = RNG_BinomialArray(Seed, Bins - 1, 0.02, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Binomial Small distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    printf("Binomial Small PMF: %.3g, %.3g, %.3g\n", RNG_BinomialPMF(0, Bins - 1, 0.02), RNG_BinomialPMF((Bins - 1) / 2, Bins - 1, 0.02), RNG_BinomialPMF(Bins - 1, Bins - 1, 0.02));

    xInt = arange(0, 1, Bins);
    FloatArray = RNG_BinomialPMFArray(xInt, Bins - 1, 0.02, Bins);
    PrintFloatArray("Binomial Small PMF distribution", FloatArray, Bins);

    free(xInt);
    free(FloatArray);

    printf("Binomial Large: " PRINT_UINT64 ", " PRINT_UINT64 ", " PRINT_UINT64 "\n", RNG_Binomial(Seed, Bins - 1, 0.98), RNG_Binomial(Seed, Bins - 1, 0.98), RNG_Binomial(NULL, Bins - 1, 0.98));

    IntArray = RNG_BinomialArray(Seed, Bins - 1, 0.98, Size);
    Hist = HistInt(0, Bins - 1, Bins, IntArray, Size);
    PrintFloatArray("Binomial Large distribution", Hist, Bins);

    free(IntArray);
    free(Hist);

    printf("Binomial Large PMF: %.3g, %.3g, %.3g\n", RNG_BinomialPMF(0, Bins - 1, 0.98), RNG_BinomialPMF((Bins - 1) / 2, Bins - 1, 0.98), RNG_BinomialPMF(Bins - 1, Bins - 1, 0.98));

    xInt = arange(0, 1, Bins);
    FloatArray = RNG_BinomialPMFArray(xInt, Bins - 1, 0.98, Bins);
    PrintFloatArray("Binomial Large PMF distribution", FloatArray, Bins);

    free(xInt);
    free(FloatArray);

    RNG_SeedDestroy(Seed);

    printf("\nErrors:\n");

    for (char *Error; (Error = RNG_ErrorArchive()) != NULL;)
        printf("%s\n", Error);

    RNG_ErrorClear();

    printf("\n");
    DBG_MemoryPrint();

    return 0;
}

void PrintIntArray(char *Name, uint64_t *Array, size_t Length)
{
    printf("%s: [", Name);

    for (uint64_t *List = Array, *EndList = Array + Length - 1; List < EndList; ++List)
        printf(PRINT_UINT64 ", ", *List);

    printf(PRINT_UINT64 "]\n", Array[Length - 1]);
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
        int64_t Index = (uint64_t)((*Numbers - Min) * Bins / (Max - Min));

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

    double Step = (Max - Min) / (double)Length;
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

double NormalPDF(double x, const void *Params)
{
    return RNG_NormalPDF(x, ((NormalParams *)Params)->mu, ((NormalParams *)Params)->sigma);
}

double ExpPDF(double x, const void *Params)
{
    return 1 / (((NormalParams *)Params)->sigma * _RNG_SQRTPI * _RNG_SQRT2) * exp(-fabs((x - ((NormalParams *)Params)->mu) / ((NormalParams *)Params)->sigma) + 0.5);
}

double ExpSampling(RNG_Seed *Seed, const void *Params)
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

double NormalSampling(RNG_Seed *Seed, double Mu, double Sigma)
{
    NormalParams Params = {.mu = Mu, .sigma = Sigma};

    return RNG_MonteCarlo(Seed, &NormalPDFArrayM, &ExpPDFArrayM, &ExpSamplingArrayM, &Params, 1);
}

double *NormalPDFArrayM(double *x, const void *Params, double *Array, size_t Size)
{
    return RNG_NormalPDFArrayM(x, ((NormalParams *)Params)->mu, ((NormalParams *)Params)->sigma, Array, Size);
}

double *ExpPDFArrayM(double *x, const void *Params, double *Array, size_t Size)
{
    double A = 1 / (((NormalParams *)Params)->sigma * _RNG_SQRTPI * _RNG_SQRT2) * exp(0.5);
    double B = 1 / ((NormalParams *)Params)->sigma;

    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-fabs(B * (*x - ((NormalParams *)Params)->mu)));

    return Array;
}

double *ExpSamplingArrayM(RNG_Seed *Seed, const void *Params, double *Array, size_t Size)
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

    return Array;
}

double *NormalSamplingArray(RNG_Seed *Seed, double Mu, double Sigma, size_t Size)
{
    NormalParams Params = {.mu = Mu, .sigma = Sigma};

    return RNG_MonteCarloArray(Seed, &NormalPDFArrayM, &ExpPDFArrayM, &ExpSamplingArrayM, &Params, 1, Size);
}
