#ifndef RANDOM2_H_INCLUDED
#define RANDOM2_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ERR_PREFIX RNG
#include <Error.h>

enum _RNG_ErrorID
{
    _RNG_ERRORID_GENERAL_NONE = 0x700000000,
    _RNG_ERRORID_INTARRAY_MALLOC = 0x700010200,
    _RNG_ERRORID_FLOATARRAY_MALLOC = 0x700020200,
    _RNG_ERRORID_EXPARRAY_MALLOC = 0x700030200,
    _RNG_ERRORID_EXPPDFARRAY_MALLOC = 0x700040200,
    _RNG_ERRORID_EXPCDFARRAY_MALLOC = 0x700050200,
    _RNG_ERRORID_EXPICDFARRAY_MALLOC = 0x700060200,
    _RNG_ERRORID_NORMALARRAY_MALLOC = 0x700030700,
    _RNG_ERRORID_NORMALPDFARRAY_MALLOC = 0x700080200,
    _RNG_ERRORID_NORMALCDFARRAY_MALLOC = 0x700090200,
    _RNG_ERRORID_NORMALICDFARRAY_MALLOC = 0x7000A0200,
    _RNG_ERRORID_GENSEED_SEED = 0x7000B0200,
    _RNG_ERRORID_CREATESEED_MALLOC = 0x7000C0200,
    _RNG_ERRORID_MONTECARLOARRAY_MALLOC = 0x7000D0200,
    _RNG_ERRORID_MONTECARLOARRAYM_MALLOC = 0x7000E0200,
    _RNG_ERRORID_MONTECARLOARRAYM_MALLOC2 = 0x7000E0201,
    _RNG_ERRORID_MONTECARLOARRAYM_MALLOC3 = 0x7000E0202
};

#define _RNG_ERRORMES_MALLOC "Unable to allocate memory (Size: %lu)"
#define _RNG_ERRORMES_CREATESEED "Unable to create seed"

typedef uint64_t* RNG_Seed;

// Generates a seed, the sed must be destroyed when no longer used
RNG_Seed RNG_GenSeed();

// Creates a seed from a number, the sed must be destroyed when no longer used
// Value: The number to convert to a seed
RNG_Seed RNG_CreateSeed(uint64_t Value);

// Destroys a seed, it must not be used after running this
void RNG_DestroySeed(RNG_Seed Seed);

// Get a uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
uint64_t RNG_Int(RNG_Seed Seed, uint64_t Min, uint64_t Max);

// Get an array of uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
uint64_t *RNG_Int_Array(RNG_Seed Seed, uint64_t Min, uint64_t Max, size_t Size);

// Fill an array with uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
void RNG_Int_ArrayM(RNG_Seed Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size);

// Get a uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
double RNG_Float(RNG_Seed Seed, double Min, double Max);

// Get an array of uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
double *RNG_Float_Array(RNG_Seed Seed, double Min, double Max, size_t Size);

// Fill an array with uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
void RNG_Float_ArrayM(RNG_Seed Seed, double Min, double Max, double *Array, size_t Size);

// Get a random number from an exponential distibution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
double RNG_Exp(RNG_Seed Seed, double x0, double l);

// Get an array of random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_Exp_Array(RNG_Seed Seed, double x0, double l, size_t Size);

// Fill an array with random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_Exp_ArrayM(RNG_Seed Seed, double x0, double l, double *Array, size_t Size);

// Get the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
double RNG_ExpPDF(double x, double x0, double l);

// Get an array of the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpPDF_Array(double *x, double x0, double l, size_t Size);

// Fill an array of the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_ExpPDF_ArrayM(double *x, double x0, double l, double *Array, size_t Size);

// Get the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
double RNG_ExpCDF(double x, double x0, double l);

// Get an array of the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpCDF_Array(double *x, double x0, double l, size_t Size);

// Fill an array of the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_ExpCDF_ArrayM(double *x, double x0, double l, double *Array, size_t Size);

// Get the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
double RNG_ExpICDF(double y, double x0, double l);

// Get an array of the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpICDF_Array(double *y, double x0, double l, size_t Size);

// Fill an array of the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_ExpICDF_ArrayM(double *y, double x0, double l, double *Array, size_t Size);

// Get a random number from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
double RNG_Normal(RNG_Seed Seed, double Mu, double Sigma);

// Get an array of random numbers from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Size: The size of the array
double *RNG_Normal_Array(RNG_Seed Seed, double Mu, double Sigma, size_t Size);

// Fill an array with random numbers from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_Normal_ArrayM(RNG_Seed Seed, double Mu, double Sigma, double *Array, size_t Size);

// Get the PDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-value to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
double RNG_NormalPDF(double x, double Mu, double Sigma);

// Get an array of the PDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Size: The size of the array
double *RNG_NormalPDF_Array(double *x, double Mu, double Sigma, size_t Size);

// Fill an array of the PDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_NormalPDF_ArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

// Get the CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-value to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
double RNG_NormalCDF(double x, double Mu, double Sigma);

// Get an array of the CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Size: The size of the array
double *RNG_NormalCDF_Array(double *x, double Mu, double Sigma, size_t Size);

// Fill an array of the CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_NormalCDF_ArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

// Get the inverse CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-value to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
double RNG_NormalICDF(double y, double Mu, double Sigma);

// Get an array of the inverse CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Size: The size of the array
double *RNG_NormalICDF_Array(double *y, double Mu, double Sigma, size_t Size);

// Fill an array of the inverse CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void *RNG_NormalICDF_ArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size);

// Poisson

// Binomial

// Samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDF: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingPDF: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <=  BoundingPDF(x), takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingSampler: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler
// Params: A pointer to some type containing all parameters for the PDFs
double RNG_MonteCarlo(RNG_Seed Seed, double (*PDF)(double x, void *Params), double (*BoundingPDF)(double x, void *Params), double (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params);

// Creates an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDF_ArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDF_ArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSampler_ArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// Size: The size of the array
double *RNG_MonteCarlo_Array(RNG_Seed Seed, void (*PDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSampler_ArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, size_t Size);

// Fills an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDF_ArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDF_ArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSampler_ArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// Size: The size of the array
// Array: The array to fill
void RNG_MonteCarlo_ArrayM(RNG_Seed Seed, void (*PDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSampler_ArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double *Array, size_t Size);

// From libit
double RNG_erfinv(double x);

// Constants
#define RNG_MAX 18446744073709551616.
#define RNG_MULTIPLIER 6364136223846793005
#define RNG_CONSTANT 1442695040888963407
#define RNG_MONTECARLO_BUFFER 1.1
#define RNG_MONTECARLO_MINSUCCESS 0.1

// For inverse error function
#define RNG_ERFINV_A3 -0.140543331
#define RNG_ERFINV_A2 0.914624893
#define RNG_ERFINV_A1 -1.645349621
#define RNG_ERFINV_A0 0.886226899

#define RNG_ERFINV_B4 0.012229801
#define RNG_ERFINV_B3 -0.329097515
#define RNG_ERFINV_B2 1.442710462
#define RNG_ERFINV_B1 -2.118377725
#define RNG_ERFINV_B0 1

#define RNG_ERFINV_C3 1.641345311
#define RNG_ERFINV_C2 3.429567803
#define RNG_ERFINV_C1 -1.62490649
#define RNG_ERFINV_C0 -1.970840454

#define RNG_ERFINV_D2 1.637067800
#define RNG_ERFINV_D1 3.543889200
#define RNG_ERFINV_D0 1

uint64_t _RNG_GlobalSeed = 0;

// Get a random uint32_t and update the seed to be that number
// Seed (RNG_Seed): The seed to use and update
#define RNG_FastInt(Seed) ((*Seed) = (*Seed) * RNG_MULTIPLIER + RNG_CONSTANT)

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (RNG_Seed): The seed to use and update
#define RNG_FastFloat(Seed) ((double) RNG_FastInt(Seed) / RNG_MAX)

RNG_Seed RNG_GenSeed()
{
    // Get the seed
    uint64_t Value = (uint64_t) time(NULL) ^ (uint64_t) clock();

    // Create seed
    RNG_Seed Seed = RNG_CreateSeed(Value);

    if (Seed == NULL)
        _RNG_AddError(_RNG_ERRORID_GENSEED_SEED, _RNG_ERRORMES_CREATESEED);

    return Seed;
}

RNG_Seed RNG_CreateSeed(uint64_t Value)
{
    // Allocate memory for it
    RNG_Seed Seed = (RNG_Seed)malloc(sizeof(uint64_t));

    if (Seed == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_CREATESEED_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(uint64_t));
        return NULL;
    }

    // Set the value
    *Seed = Value;

    return Seed;
}

void RNG_DestroySeed(RNG_Seed Seed)
{
    free(Seed);
}

uint64_t RNG_Int(uint64_t *Seed, uint64_t Min, uint64_t Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + (RNG_FastInt(Seed) % (1 + Max - Min));
}

uint64_t *RNG_Int_Array(uint64_t *Seed, uint64_t Min, uint64_t Max, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_INTARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    RNG_Int_ArrayM(Seed, Min, Max, Array, Size);

    return Array;
}

void RNG_Int_ArrayM(uint64_t *Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + (RNG_FastInt(Seed) % (1 + Max - Min));
}

double RNG_Float(uint64_t *Seed, double Min, double Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + RNG_FastFloat(Seed) * (Max - Min);
}

double *RNG_Float_Array(uint64_t *Seed, double Min, double Max, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_FLOATARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_Float_ArrayM(Seed, Min, Max, Array, Size);

    return Array;
}

void RNG_Float_ArrayM(uint64_t *Seed, double Min, double Max, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + RNG_FastFloat(Seed) * (Max - Min);
}

double RNG_Exp(uint64_t *Seed, double x0, double l)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FastFloat(Seed);

    // Get from current distribution
    return x0 - l * log(1 - Uniform);
}

double *RNG_Exp_Array(uint64_t *Seed, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_Exp_ArrayM(Seed, x0, l, Array, Size);

    return Array;
}

void RNG_Exp_ArrayM(uint64_t *Seed, double x0, double l, double *Array, size_t Size)
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

        // Get from distribution
        *List = x0 - l * log(1 - Uniform);
    }
}

double RNG_ExpPDF(double x, double x0, double l)
{
    // Calculate constants
    double c = 1 / l;

    // Calculate the PDF
    return ((x < x0) ? (0) : (c * exp(-c * (x - x0))));
}

double *RNG_ExpPDF_Array(double *x, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPPDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpPDF_ArrayM(x, x0, l, Array, Size);

    return Array;
}

void *RNG_ExpPDF_ArrayM(double *x, double x0, double l, double *Array, size_t Size)
{
    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < x0) ? (0) : (c * exp(-c * (*x - x0))));
}

double RNG_ExpCDF(double x, double x0, double l)
{
    // Calculate constants
    double c = 1 / l;

    // Calculate the CDF
    return ((x < x0) ? (0) : (1 - exp(-c * (x - x0))));
}

double *RNG_ExpCDF_Array(double *x, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPCDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpCDF_ArrayM(x, x0, l, Array, Size);

    return Array;
}

void *RNG_ExpCDF_ArrayM(double *x, double x0, double l, double *Array, size_t Size)
{
    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < x0) ? (0) : (1 - exp(-c * (*x - x0))));
}

double RNG_ExpICDF(double y, double x0, double l)
{
    // Calculate the ICDF
    return x0 - l * log(1 - y);
}

double *RNG_ExpICDF_Array(double *y, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPICDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpICDF_ArrayM(y, x0, l, Array, Size);

    return Array;
}

void *RNG_ExpICDF_ArrayM(double *y, double x0, double l, double *Array, size_t Size)
{
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = x0 - l * log(1 - *y);
}

double RNG_Normal(uint64_t *Seed, double Mu, double Sigma)
{
    // Calulate constants
    double A = M_SQRT2 * Sigma;

    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FastFloat(Seed);

    // Get from current distribution
    return Mu + A * RNG_erfinv(2 * Uniform - 1);
}

double *RNG_Normal_Array(uint64_t *Seed, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_NORMALARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_Normal_ArrayM(Seed, Mu, Sigma, Array, Size);

    return Array;
}

void RNG_Normal_ArrayM(uint64_t *Seed, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calulate constants
    double A = M_SQRT2 * Sigma;

    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FastFloat(Seed);

        // Get from distribution
        *List = Mu + A * RNG_erfinv(2 * Uniform - 1);
    }
}

double RNG_NormalPDF(double x, double Mu, double Sigma)
{
    // Calculate constants
    double A = M_2_SQRTPI / (2 * Sigma * M_SQRT2);
    double B = 0.5 / (Sigma * Sigma);

    // Calculate the PDF
    return A * exp(-B * (x - Mu) * (x - Mu));
}

double *RNG_NormalPDF_Array(double *x, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPPDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalPDF_ArrayM(x, Mu, Sigma, Array, Size);

    return Array;
}

void *RNG_NormalPDF_ArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = M_2_SQRTPI / (2 * Sigma * M_SQRT2);
    double B = 0.5 / (Sigma * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-B * (*x - Mu) * (*x - Mu));
}

double RNG_NormalCDF(double x, double Mu, double Sigma)
{
    // Calculate constants
    double A = 1 / (M_SQRT2 * Sigma);

    // Calculate the CDF
    return 0.5 * (1 + erf(A * (x - Mu)));
}

double *RNG_NormalCDF_Array(double *x, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPCDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalCDF_ArrayM(x, Mu, Sigma, Array, Size);

    return Array;
}

void *RNG_NormalCDF_ArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = 1 / (M_SQRT2 * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = 0.5 * (1 + erf(A * (*x - Mu)));
}

double RNG_NormalICDF(double y, double Mu, double Sigma)
{
    // Calculate constants
    double A = 2 * Sigma;

    // Calculate the ICDF
    return Mu + A * RNG_erfinv(2 * y - 1);
}

double *RNG_NormalICDF_Array(double *y, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPICDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalICDF_ArrayM(y, Mu, Sigma, Array, Size);

    return Array;
}

void *RNG_NormalICDF_ArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = 2 * Sigma;
 
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = Mu + A * RNG_erfinv(2 * *y - 1);
}

double RNG_MonteCarlo(RNG_Seed Seed, double (*PDF)(double x, void *Params), double (*BoundingPDF)(double x, void *Params), double (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Loop until a value has been found
    while (true)
    {
        // Get a random sample
        double x = (*BoundingSampler)(Seed, Params);

        // Get the success probability
        double Prob = (*PDF)(x, Params) / (*BoundingPDF)(x, Params);

        // Check to accept it
        double Uniform = RNG_FastFloat(Seed);

        if (Uniform < Prob)
            return x;
    }
}

double *RNG_MonteCarlo_Array(RNG_Seed Seed, void (*PDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSampler_ArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_MONTECARLOARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_MonteCarlo_ArrayM(Seed, PDF_ArrayM, BoundingPDF_ArrayM, BoundingSampler_ArrayM, Params, Array, Size);

    return Array;
}

void RNG_MonteCarlo_ArrayM(RNG_Seed Seed, void (*PDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDF_ArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSampler_ArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    double *x = (double *)malloc(sizeof(double) * Size);

    if (x == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_MONTECARLOARRAYM_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return;
    }

    double *PDF = (double *)malloc(sizeof(double) * Size);

    if (PDF == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_MONTECARLOARRAYM_MALLOC2, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(x);
        return;
    }

    double *BoundingPDF = (double *)malloc(sizeof(double) * Size);

    if (BoundingPDF == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_MONTECARLOARRAYM_MALLOC3, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(x);
        free(PDF);
        return;
    }

    // Set end point
    double *ArrayEnd = Array + Size;
    double SuccessRate = 1;

    // Loop until done
    while (true)
    {
        // Get x samples
        size_t SampleCount = (size_t)((double)(ArrayEnd - Array) / SuccessRate * RNG_MONTECARLO_BUFFER);

        if (SampleCount > Size)
            SampleCount = Size;

        BoundingSampler_ArrayM(Seed, Params, x, SampleCount);

        // Get the PDF values
        PDF_ArrayM(x, Params, PDF, SampleCount);
        BoundingPDF_ArrayM(x, Params, BoundingPDF, SampleCount);

        // Check if they are good
        double *OldArray = Array;

        for (double *xSample = x, *PDFSample = PDF, *BoundingPDFSample = BoundingPDF, *xSampleEnd = xSample + SampleCount; xSample < xSampleEnd; ++xSample, ++PDFSample, ++BoundingPDFSample)
        {
            // Get probability for success
            double Prob = *PDFSample / *BoundingPDFSample;
            double Uniform = RNG_FastFloat(Seed);

            // Check if it should be kept
            if (Uniform < Prob)
            {
                *(Array++) = *xSample;

                // Check if it is done
                if (Array >= ArrayEnd)
                {
                    free(x);
                    free(PDF);
                    free(BoundingPDF);

                    return;
                }
            }
        }

        // Get the new success rate
        SuccessRate = (double)(Array - OldArray) / (double)SampleCount;

        if (SuccessRate < RNG_MONTECARLO_MINSUCCESS)
            SuccessRate = RNG_MONTECARLO_MINSUCCESS;
    }
}

// From libit
double RNG_erfinv(double x)
{
    double x2, r, y;
    int32_t sign_x;

    if (x < -1 || x > 1)
        return NAN;

    if (x == 0)
        return 0;

    if (x > 0)
        sign_x = 1;

    else 
    {
        sign_x = -1;
        x = -x;
    }

    if (x <= 0.7) 
    {
        x2 = x * x;
        r = x * (((RNG_ERFINV_A3 * x2 + RNG_ERFINV_A2) * x2 + RNG_ERFINV_A1) * x2 + RNG_ERFINV_A0);
        r /= (((RNG_ERFINV_B4 * x2 + RNG_ERFINV_B3) * x2 + RNG_ERFINV_B2) * x2 + RNG_ERFINV_B1) * x2 + RNG_ERFINV_B0;
    }
  
    else 
    {
        y = sqrt (-log ((1 - x) / 2));
        r = (((RNG_ERFINV_C3 * y + RNG_ERFINV_C2) * y + RNG_ERFINV_C1) * y + RNG_ERFINV_C0);
        r /= ((RNG_ERFINV_D2 * y + RNG_ERFINV_D1) * y + RNG_ERFINV_D0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf(r) - x) / (M_2_SQRTPI * exp(-r * r));
    r -= (erf(r) - x) / (M_2_SQRTPI * exp(-r * r));

    return r;
}

#endif // RANDOM2_H_INCLUDED
