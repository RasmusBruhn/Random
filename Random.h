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
    _RNG_ERRORID_EXPICDFARRAY_MALLOC = 0x700060200
};

#define _RNG_ERRORMES_MALLOC "Unable to allocate memory (Size: %lu)"

// Get a uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
uint64_t RNG_Int(uint64_t *Seed, uint64_t Min, uint64_t Max);

// Get an array of uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
uint64_t *RNG_Int_Array(uint64_t *Seed, uint64_t Min, uint64_t Max, size_t Size);

// Fill an array with uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array The array to fill
// Size: The size of the array
void RNG_Int_ArrayM(uint64_t *Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size);

// Get a uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
double RNG_Float(uint64_t *Seed, double Min, double Max);

// Get an array of uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
double *RNG_Float_Array(uint64_t *Seed, double Min, double Max, size_t Size);

// Fill an array with uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array The array to fill
// Size: The size of the array
void RNG_Float_ArrayM(uint64_t *Seed, double Min, double Max, double *Array, size_t Size);

// Get a random number from an exponential distibution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// l: The length scale of the distribution
double RNG_Exp(uint64_t *Seed, double l);

// Get an array of random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_Exp_Array(uint64_t *Seed, double l, size_t Size);

// Fill an array with random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// l: The length scale of the distribution
// Array The array to fill
// Size: The size of the array
void RNG_Exp_ArrayM(uint64_t *Seed, double l, double *Array, size_t Size);

// Get the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// l: The length scale of the distribution
double RNG_ExpPDF(double x, double l);

// Get an array of the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpPDF_Array(double *x, double l, size_t Size);

// Fill an array of the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Array The array to fill
// Size: The size of the array
void *RNG_ExpPDF_ArrayM(double *x, double l, double *Array, size_t Size);

// Get the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// l: The length scale of the distribution
double RNG_ExpCDF(double x, double l);

// Get an array of the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpCDF_Array(double *x, double l, size_t Size);

// Fill an array of the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Array The array to fill
// Size: The size of the array
void *RNG_ExpCDF_ArrayM(double *x, double l, double *Array, size_t Size);

// Get the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-value to get the y-value for
// l: The length scale of the distribution
double RNG_ExpICDF(double y, double l);

// Get an array of the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpICDF_Array(double *y, double l, size_t Size);

// Fill an array of the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// l: The length scale of the distribution
// Array The array to fill
// Size: The size of the array
void *RNG_ExpICDF_ArrayM(double *y, double l, double *Array, size_t Size);

// DoubleExp

// Normal

// Poisson

// Binomial

// Monte carlo

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

uint64_t RNG_Int(uint64_t *Seed, uint64_t Min, uint64_t Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + (RNG_FastInt(*Seed) % (1 + Max - Min));
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
        *List = Min + (RNG_FastInt(*Seed) % (1 + Max - Min));
}

double RNG_Float(uint64_t *Seed, double Min, double Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + RNG_FastFloat(*Seed) * (Max - Min);
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
        *List = Min + RNG_FastFloat(*Seed) * (Max - Min);
}

double RNG_Exp(uint64_t *Seed, double l)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FastFloat(*Seed);

    // Get from current distribution
    return -l * log(1 - Uniform);
}

double *RNG_Exp_Array(uint64_t *Seed, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_Exp_ArrayM(Seed, l, Array, Size);

    return Array;
}

void RNG_Exp_ArrayM(uint64_t *Seed, double l, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FastFloat(*Seed);

        // Get from distribution
        *List = -l * log(1 - Uniform);
    }
}

double RNG_ExpPDF(double x, double l)
{
    // Calculate constants
    double c = 1 / l;

    // Calculate the PDF
    return ((x < 0) ? (0) : (c * exp(-c * x)));
}

double *RNG_ExpPDF_Array(double *x, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPPDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpPDF_ArrayM(x, l, Array, Size);

    return Array;
}

void *RNG_ExpPDF_ArrayM(double *x, double l, double *Array, size_t Size)
{
    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < 0) ? (0) : (c * exp(-c * *x)));
}

double RNG_ExpCDF(double x, double l)
{
    // Calculate constants
    double c = 1 / l;

    // Calculate the CDF
    return ((x < 0) ? (0) : (1 - exp(-c * x)));
}

double *RNG_ExpCDF_Array(double *x, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPCDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpCDF_ArrayM(x, l, Array, Size);

    return Array;
}

void *RNG_ExpCDF_ArrayM(double *x, double l, double *Array, size_t Size)
{
    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < 0) ? (0) : (1 - exp(-c * *x)));
}

double RNG_ExpICDF(double y, double l)
{
    // Calculate the ICDF
    return -l * log(1 - y);
}

double *RNG_ExpICDF_Array(double *y, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_SetError(_RNG_ERRORID_EXPICDFARRAY_MALLOC, _RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpICDF_ArrayM(y, l, Array, Size);

    return Array;
}

void *RNG_ExpICDF_ArrayM(double *y, double l, double *Array, size_t Size)
{
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = -l * log(1 - *y);
}

#endif // RANDOM2_H_INCLUDED
