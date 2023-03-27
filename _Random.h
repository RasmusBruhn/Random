#ifndef RANDOM_H_INCLUDED2
#define RANDOM_H_INCLUDED2

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Random.h"

#define ERR_PREFIX RNG
#include <Error2.c>
#define _RNG_ErrorSet(Format, ...) __RNG_ErrorSet(__FILE__, __LINE__, Format __VA_OPT__(, ) __VA_ARGS__)
#define _RNG_ErrorAdd(Format, ...) __RNG_ErrorAdd(__FILE__, __LINE__, Format __VA_OPT__(, ) __VA_ARGS__)
#define _RNG_ErrorAddExternal(ExternalMessage, Format, ...) __RNG_ErrorAddExternal(__FILE__, __LINE__, ExternalMessage, Format __VA_OPT__(, ) __VA_ARGS__)

#define _RNG_ERRORMES_MALLOC "Unable to allocate memory (Size: %lu)"
#define _RNG_ERRORMES_CREATESEED "Unable to create seed"
#define _RNG_ERRORMES_GENERATE "Unable to generate random %s"
#define _RNG_ERRORMES_PMF "Unable to generate PMF for %s"
#define _RNG_ERRORMES_PDF "Unable to generate PDF for %s"
#define _RNG_ERRORMES_CDF "Unable to generate CDF for %s"
#define _RNG_ERRORMES_ICDF "Unable to generate ICDF for %s"
#define _RNG_ERRORMES_LARGERTHAN "The value of %s must not be smaller than %s but received (%lu, %lu)"
#define _RNG_ERRORMES_LARGERTHANF "The value of %s must not be smaller than %s but received (%.3g, %.3g)"
#define _RNG_ERRORMES_NOTNEGATIVE "The value of %s must not be negative, but received %.3g"
#define _RNG_ERRORMES_INRANGE "The value of %s must be in the range [%.3g, %.3g] but received %.3g"
#define _RNG_ERRORMES_MONTECARLOSAMPLE "Unable to sample from bounding sampler"
#define _RNG_ERRORMES_MONTECARLOPDF "Unable to evaluate PDF value"
#define _RNG_ERRORMES_MONTECARLOBOUNDINGPDF "Unable to evaluate bounding PDF value"
#define _RNG_ERRORMES_MONTECARLOPMF "Unable to evaluate PMF value"
#define _RNG_ERRORMES_MONTECARLOBOUNDINGPMF "Unable to evaluate bounding PMF value"

typedef struct ___RNG_BinomialParams _RNG_BinomialParams;

struct ___RNG_BinomialParams
{
    uint64_t N;
    double p;
};

// Fills an array of PMF values for the poisson distribution using params (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_PoissonPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size);

// Fills an array of random numbers from the poisson bounding distribution
// Seed: The seed to use and update, NULL to use global seed
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
uint64_t *_RNG_PoissonBoundingSamplerArrayM(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size);

// Fills an array of bounding PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_PoissonBoundingPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size);

// Fills an array of PMF values for the binomial distribution using params (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_BinomialPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size);

// Fills an array of random numbers from the binomial bounding distribution
// Seed: The seed to use and update, NULL to use global seed
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
uint64_t *_RNG_BinomialBoundingSamplerArrayM(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size);

// Fills an array of bounding PMF values for the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_BinomialBoundingPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size);

// From libit
double _RNG_erfinv(double x);

// Constants
#define _RNG_MONTECARLO_BUFFER 1.1
#define _RNG_MONTECARLO_MINSUCCESS 0.1

// For inverse error function
#define _RNG_ERFINV_A3 -0.140543331
#define _RNG_ERFINV_A2 0.914624893
#define _RNG_ERFINV_A1 -1.645349621
#define _RNG_ERFINV_A0 0.886226899

#define _RNG_ERFINV_B4 0.012229801
#define _RNG_ERFINV_B3 -0.329097515
#define _RNG_ERFINV_B2 1.442710462
#define _RNG_ERFINV_B1 -2.118377725
#define _RNG_ERFINV_B0 1

#define _RNG_ERFINV_C3 1.641345311
#define _RNG_ERFINV_C2 3.429567803
#define _RNG_ERFINV_C1 -1.62490649
#define _RNG_ERFINV_C0 -1.970840454

#define _RNG_ERFINV_D2 1.637067800
#define _RNG_ERFINV_D1 3.543889200
#define _RNG_ERFINV_D0 1

// Constants
#define _RNG_SQRT2 1.4142135623730951
#define _RNG_SQRTPI 1.7724538509055159
#define _RNG_1_E 0.36787944117144233
#define _RNG_E 2.718281828459045

#endif