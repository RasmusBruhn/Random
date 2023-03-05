#ifndef RANDOM2_H_INCLUDED
#define RANDOM2_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ERR_PREFIX RNG
#include <Error2.h>
#define _RNG_ErrorSet(Format, ...) __RNG_ErrorSet(__FILE__, __LINE__, Format __VA_OPT__(, ) __VA_ARGS__)
#define _RNG_ErrorAdd(Format, ...) __RNG_ErrorAdd(__FILE__, __LINE__, Format __VA_OPT__(, ) __VA_ARGS__)
#define _RNG_ErrorAddExternal(ExternalMessage, Format, ...) __RNG_ErrorAddExternal(__FILE__, __LINE__, ExternalMessage, Format __VA_OPT__(, ) __VA_ARGS__)

#define _RNG_ERRORMES_MALLOC "Unable to allocate memory (Size: %lu)"
#define _RNG_ERRORMES_CREATESEED "Unable to create seed"

typedef uint64_t* RNG_Seed;

// Generates a seed, the sed must be destroyed when no longer used
RNG_Seed RNG_SeedGenerate();

// Creates a seed from a number, the sed must be destroyed when no longer used
// Value: The number to convert to a seed
RNG_Seed RNG_SeedCreate(uint64_t Value);

// Destroys a seed, it must not be used after running this
void RNG_SeedDestroy(RNG_Seed Seed);

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
uint64_t *RNG_IntArray(RNG_Seed Seed, uint64_t Min, uint64_t Max, size_t Size);

// Fill an array with uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
void RNG_IntArrayM(RNG_Seed Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size);

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
double *RNG_FloatArray(RNG_Seed Seed, double Min, double Max, size_t Size);

// Fill an array with uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
void RNG_FloatArrayM(RNG_Seed Seed, double Min, double Max, double *Array, size_t Size);

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
double *RNG_ExpArray(RNG_Seed Seed, double x0, double l, size_t Size);

// Fill an array with random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_ExpArrayM(RNG_Seed Seed, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpPDFArray(double *x, double x0, double l, size_t Size);

// Fill an array of the PDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_ExpPDFArrayM(double *x, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpCDFArray(double *x, double x0, double l, size_t Size);

// Fill an array of the CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_ExpCDFArrayM(double *x, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpICDFArray(double *y, double x0, double l, size_t Size);

// Fill an array of the inverse CDF for an exponential distribution (1 / l * exp(-x / l))
// x: The x-values to get the y-value for
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_ExpICDFArrayM(double *y, double x0, double l, double *Array, size_t Size);

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
double *RNG_NormalArray(RNG_Seed Seed, double Mu, double Sigma, size_t Size);

// Fill an array with random numbers from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_NormalArrayM(RNG_Seed Seed, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalPDFArray(double *x, double Mu, double Sigma, size_t Size);

// Fill an array of the PDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_NormalPDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalCDFArray(double *x, double Mu, double Sigma, size_t Size);

// Fill an array of the CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_NormalCDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalICDFArray(double *y, double Mu, double Sigma, size_t Size);

// Fill an array of the inverse CDF for a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// x: The x-values to get the y-value for
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_NormalICDFArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size);

// Get a random number from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
uint32_t RNG_Poisson(RNG_Seed Seed, double Lambda);

// Get an array of random numbers from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
// Size: The size of the array
uint32_t *RNG_PoissonArray(RNG_Seed Seed, double Lambda, size_t Size);

// Fills an array of random numbers from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
// Array: The array to fill
// Size: The size of the array
void RNG_PoissonArrayM(RNG_Seed Seed, double Lambda, uint32_t *Array, size_t Size);

// Get the PMF for a poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The position to get the PMF for
// Lambda: The mean value of the distribution
double RNG_PoissonPMF(uint32_t n, double Lambda);

// Get an array of PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Lambda: The mean value of the distribution
// Size: The size of the array
double *RNG_PoissonPMFArray(uint32_t *n, double Lambda, size_t Size);

// Fills an array of PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Lambda: The mean value of the distribution
// Array: The array to fill
// Size: The size of the array
double *RNG_PoissonPMFArrayM(uint32_t *n, double Lambda, double *Array, size_t Size);

// Get the PMF for a poisson distribution using params (lambda ** n * exp(lambda) / n!)
// n: The position to get the PMF for
// Params: The parameters
double _RNG_PoissonPMF(uint32_t n, void *Params);

// Get an array of PMF values for the poisson distribution using params (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Size: The size of the array
double *_RNG_PoissonPMFArray(uint32_t *n, void *Params, size_t Size);

// Fills an array of PMF values for the poisson distribution using params (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_PoissonPMFArrayM(uint32_t *n, void *Params, double *Array, size_t Size);

// Get a random number from the poisson bounding distribution
// Seed: The seed to use and update, NULL to use global seed
// Params: The parameters
uint32_t _RNG_PoissonBoundingSampler(RNG_Seed Seed, void *Params);

// Get an array of random numbers from the poisson bounding distribution
// Seed: The seed to use and update, NULL to use global seed
// Params: The parameters
// Size: The size of the array
uint32_t *_RNG_PoissonBoundingSamplerArray(RNG_Seed Seed, void *Params, size_t Size);

// Fills an array of random numbers from the poisson bounding distribution
// Seed: The seed to use and update, NULL to use global seed
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
void _RNG_PoissonBoundingSamplerArrayM(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size);

// Get the bounding PMF for a poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The position to get the PMF for
// Params: The parameters
double _RNG_PoissonBoundingPMF(uint32_t n, void *Params);

// Get an array of bounding PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Size: The size of the array
double *_RNG_PoissonBoundingPMFArray(uint32_t *n, void *Params, size_t Size);

// Fills an array of bounding PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Params: The parameters
// Array: The array to fill
// Size: The size of the array
double *_RNG_PoissonBoundingPMFArrayM(uint32_t *n, void *Params, double *Array, size_t Size);

// Binomial

// Samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDF: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingPDF: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingSampler: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
double RNG_MonteCarlo(RNG_Seed Seed, double (*PDF)(double x, void *Params), double (*BoundingPDF)(double x, void *Params), double (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params, double BoundingMultiplier);

// Creates an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
double *RNG_MonteCarloArray(RNG_Seed Seed, void (*PDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double BoundingMultiplier, size_t Size);

// Fills an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
// Array: The array to fill
void RNG_MonteCarloArrayM(RNG_Seed Seed, void (*PDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double BoundingMultiplier, double *Array, size_t Size);

// Samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMF: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingPMF: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingSampler: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
uint32_t RNG_MonteCarloUInt(RNG_Seed Seed, double (*PMF)(uint32_t n, void *Params), double (*BoundingPMF)(uint32_t n, void *Params), uint32_t (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params, double BoundingMultiplier);

// Creates an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPMFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
uint32_t *RNG_MonteCarloUIntArray(RNG_Seed Seed, void (*PMFArrayM)(uint32_t *n, void *Params, double *Array, size_t Size), void (*BoundingPMFArrayM)(uint32_t *n, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size), void *Params, double BoundingMultiplier, size_t Size);

// Fills an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPMFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
// Array: The array to fill
void RNG_MonteCarloUIntArrayM(RNG_Seed Seed, void (*PMFArrayM)(uint32_t *x, void *Params, double *Array, size_t Size), void (*BoundingPMFArrayM)(uint32_t *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size), void *Params, double BoundingMultiplier, uint32_t *Array, size_t Size);

// From libit
double _RNG_erfinv(double x);

// Constants
#define RNG_MAX 18446744073709551616.
#define _RNG_MULTIPLIER 6364136223846793005
#define _RNG_CONSTANT 1442695040888963407
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

uint64_t _RNG_GlobalSeed = 0;

// Get a random uint32_t and update the seed to be that number
// Seed (RNG_Seed): The seed to use and update
#define RNG_IntFast(Seed) ((*Seed) = (*Seed) * _RNG_MULTIPLIER + _RNG_CONSTANT)

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (RNG_Seed): The seed to use and update
#define RNG_FloatFast(Seed) ((double) RNG_IntFast(Seed) / RNG_MAX)

RNG_Seed RNG_SeedGenerate()
{
    // Get the seed
    uint64_t Value = (uint64_t) time(NULL) ^ (uint64_t) clock();

    // Create seed
    RNG_Seed Seed = RNG_SeedCreate(Value);

    if (Seed == NULL)
        _RNG_ErrorAdd(_RNG_ERRORMES_CREATESEED);

    return Seed;
}

RNG_Seed RNG_SeedCreate(uint64_t Value)
{
    // Allocate memory for it
    RNG_Seed Seed = (RNG_Seed)malloc(sizeof(uint64_t));

    if (Seed == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t));
        return NULL;
    }

    // Set the value
    *Seed = Value;

    return Seed;
}

void RNG_SeedDestroy(RNG_Seed Seed)
{
    free(Seed);
}

uint64_t RNG_Int(RNG_Seed Seed, uint64_t Min, uint64_t Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + (RNG_IntFast(Seed) % (1 + Max - Min));
}

uint64_t *RNG_IntArray(RNG_Seed Seed, uint64_t Min, uint64_t Max, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    RNG_IntArrayM(Seed, Min, Max, Array, Size);

    return Array;
}

void RNG_IntArrayM(RNG_Seed Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + (RNG_IntFast(Seed) % (1 + Max - Min));
}

double RNG_Float(RNG_Seed Seed, double Min, double Max)
{
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    return Min + RNG_FloatFast(Seed) * (Max - Min);
}

double *RNG_FloatArray(RNG_Seed Seed, double Min, double Max, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_FloatArrayM(Seed, Min, Max, Array, Size);

    return Array;
}

void RNG_FloatArrayM(RNG_Seed Seed, double Min, double Max, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + RNG_FloatFast(Seed) * (Max - Min);
}

double RNG_Exp(RNG_Seed Seed, double x0, double l)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FloatFast(Seed);

    // Get from current distribution
    return x0 - l * log(1 - Uniform);
}

double *RNG_ExpArray(RNG_Seed Seed, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpArrayM(Seed, x0, l, Array, Size);

    return Array;
}

void RNG_ExpArrayM(RNG_Seed Seed, double x0, double l, double *Array, size_t Size)
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

double *RNG_ExpPDFArray(double *x, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpPDFArrayM(x, x0, l, Array, Size);

    return Array;
}

void RNG_ExpPDFArrayM(double *x, double x0, double l, double *Array, size_t Size)
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

double *RNG_ExpCDFArray(double *x, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpCDFArrayM(x, x0, l, Array, Size);

    return Array;
}

void RNG_ExpCDFArrayM(double *x, double x0, double l, double *Array, size_t Size)
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

double *RNG_ExpICDFArray(double *y, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_ExpICDFArrayM(y, x0, l, Array, Size);

    return Array;
}

void RNG_ExpICDFArrayM(double *y, double x0, double l, double *Array, size_t Size)
{
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = x0 - l * log(1 - *y);
}

double RNG_Normal(RNG_Seed Seed, double Mu, double Sigma)
{
    // Calulate constants
    double A = _RNG_SQRT2 * Sigma;

    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get the number from the uniform distribution
    double Uniform = RNG_FloatFast(Seed);

    // Get from current distribution
    return Mu + A * _RNG_erfinv(2 * Uniform - 1);
}

double *RNG_NormalArray(RNG_Seed Seed, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalArrayM(Seed, Mu, Sigma, Array, Size);

    return Array;
}

void RNG_NormalArrayM(RNG_Seed Seed, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calulate constants
    double A = _RNG_SQRT2 * Sigma;

    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FloatFast(Seed);

        // Get from distribution
        *List = Mu + A * _RNG_erfinv(2 * Uniform - 1);
    }
}

double RNG_NormalPDF(double x, double Mu, double Sigma)
{
    // Calculate constants
    double A = 1 / (Sigma * _RNG_SQRTPI * _RNG_SQRT2);
    double B = 0.5 / (Sigma * Sigma);

    // Calculate the PDF
    return A * exp(-B * (x - Mu) * (x - Mu));
}

double *RNG_NormalPDFArray(double *x, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalPDFArrayM(x, Mu, Sigma, Array, Size);

    return Array;
}

void RNG_NormalPDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = 1 / (Sigma * _RNG_SQRTPI * _RNG_SQRT2);
    double B = 0.5 / (Sigma * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-B * (*x - Mu) * (*x - Mu));
}

double RNG_NormalCDF(double x, double Mu, double Sigma)
{
    // Calculate constants
    double A = 1 / (_RNG_SQRT2 * Sigma);

    // Calculate the CDF
    return 0.5 * (1 + erf(A * (x - Mu)));
}

double *RNG_NormalCDFArray(double *x, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalCDFArrayM(x, Mu, Sigma, Array, Size);

    return Array;
}

void RNG_NormalCDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = 1 / (_RNG_SQRT2 * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = 0.5 * (1 + erf(A * (*x - Mu)));
}

double RNG_NormalICDF(double y, double Mu, double Sigma)
{
    // Calculate constants
    double A = 2 * Sigma;

    // Calculate the ICDF
    return Mu + A * _RNG_erfinv(2 * y - 1);
}

double *RNG_NormalICDFArray(double *y, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_NormalICDFArrayM(y, Mu, Sigma, Array, Size);

    return Array;
}

void RNG_NormalICDFArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size)
{
    // Calculate constants
    double A = 2 * Sigma;
 
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = Mu + A * _RNG_erfinv(2 * *y - 1);
}

double RNG_MonteCarlo(RNG_Seed Seed, double (*PDF)(double x, void *Params), double (*BoundingPDF)(double x, void *Params), double (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params, double BoundingMultiplier)
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
        double Prob = (*PDF)(x, Params) / (BoundingMultiplier * (*BoundingPDF)(x, Params));

        // Check to accept it
        double Uniform = RNG_FloatFast(Seed);

        if (Uniform < Prob)
            return x;
    }
}

double *RNG_MonteCarloArray(RNG_Seed Seed, void (*PDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double BoundingMultiplier, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_MonteCarloArrayM(Seed, PDFArrayM, BoundingPDFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, Array, Size);

    return Array;
}

void RNG_MonteCarloArrayM(RNG_Seed Seed, void (*PDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingPDFArrayM)(double *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, double *Array, size_t Size), void *Params, double BoundingMultiplier, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    double *x = (double *)malloc(sizeof(double) * Size);

    if (x == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return;
    }

    double *PDF = (double *)malloc(sizeof(double) * Size);

    if (PDF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(x);
        return;
    }

    double *BoundingPDF = (double *)malloc(sizeof(double) * Size);

    if (BoundingPDF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
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
        size_t SampleCount = (size_t)((double)(ArrayEnd - Array) / SuccessRate * _RNG_MONTECARLO_BUFFER);

        if (SampleCount > Size)
            SampleCount = Size;

        BoundingSamplerArrayM(Seed, Params, x, SampleCount);

        // Get the PDF values
        PDFArrayM(x, Params, PDF, SampleCount);
        BoundingPDFArrayM(x, Params, BoundingPDF, SampleCount);

        // Check if they are good
        double *OldArray = Array;

        for (double *xSample = x, *PDFSample = PDF, *BoundingPDFSample = BoundingPDF, *xSampleEnd = xSample + SampleCount; xSample < xSampleEnd; ++xSample, ++PDFSample, ++BoundingPDFSample)
        {
            // Get probability for success
            double Prob = *PDFSample / (BoundingMultiplier * *BoundingPDFSample);
            double Uniform = RNG_FloatFast(Seed);

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

        if (SuccessRate < _RNG_MONTECARLO_MINSUCCESS)
            SuccessRate = _RNG_MONTECARLO_MINSUCCESS;
    }
}

uint32_t RNG_MonteCarloUInt(RNG_Seed Seed, double (*PMF)(uint32_t n, void *Params), double (*BoundingPMF)(uint32_t n, void *Params), uint32_t (*BoundingSampler)(RNG_Seed Seed, void *Params), void *Params, double BoundingMultiplier)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Loop until a value has been found
    while (true)
    {
        // Get a random sample
        uint32_t n = (uint32_t)(*BoundingSampler)(Seed, Params);

        // Get the success probability
        double Prob = (*PMF)(n, Params) / (BoundingMultiplier * (*BoundingPMF)(n, Params));

        // Check to accept it
        double Uniform = RNG_FloatFast(Seed);

        if (Uniform < Prob)
            return n;
    }
}

uint32_t *RNG_MonteCarloUIntArray(RNG_Seed Seed, void (*PMFArrayM)(uint32_t *n, void *Params, double *Array, size_t Size), void (*BoundingPMFArrayM)(uint32_t *n, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size), void *Params, double BoundingMultiplier, size_t Size)
{
    // Get memory
    uint32_t *Array = (uint32_t *)malloc(sizeof(uint32_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint32_t) * Size);
        return NULL;
    }

    // Get numbers
    RNG_MonteCarloUIntArrayM(Seed, PMFArrayM, BoundingPMFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, Array, Size);

    return Array;
}

void RNG_MonteCarloUIntArrayM(RNG_Seed Seed, void (*PMFArrayM)(uint32_t *n, void *Params, double *Array, size_t Size), void (*BoundingPMFArrayM)(uint32_t *x, void *Params, double *Array, size_t Size), void (*BoundingSamplerArrayM)(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size), void *Params, double BoundingMultiplier, uint32_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    uint32_t *n = (uint32_t *)malloc(sizeof(uint32_t) * Size);

    if (n == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint32_t) * Size);
        return;
    }

    double *PMF = (double *)malloc(sizeof(double) * Size);

    if (PMF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(n);
        return;
    }

    double *BoundingPMF = (double *)malloc(sizeof(double) * Size);

    if (BoundingPMF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(n);
        free(PMF);
        return;
    }

    // Set end point
    double *ArrayEnd = Array + Size;
    double SuccessRate = 1;

    // Loop until done
    while (true)
    {
        // Get x samples
        size_t SampleCount = (size_t)((double)(ArrayEnd - Array) / SuccessRate * _RNG_MONTECARLO_BUFFER);

        if (SampleCount > Size)
            SampleCount = Size;

        BoundingSamplerArrayM(Seed, Params, n, SampleCount);

        // Get the PDF values
        PMFArrayM(n, Params, PMF, SampleCount);
        BoundingPMFArrayM(n, Params, BoundingPMF, SampleCount);

        // Check if they are good
        double *OldArray = Array;

        for (double *nSample = n, *PMFSample = PMF, *BoundingPMFSample = BoundingPMF, *nSampleEnd = nSample + SampleCount; nSample < nSampleEnd; ++nSample, ++PMFSample, ++BoundingPMFSample)
        {
            // Get probability for success
            double Prob = *PMFSample / (BoundingMultiplier * *BoundingPMFSample);
            double Uniform = RNG_FloatFast(Seed);

            // Check if it should be kept
            if (Uniform < Prob)
            {
                *(Array++) = *nSample;

                // Check if it is done
                if (Array >= ArrayEnd)
                {
                    free(n);
                    free(PMF);
                    free(BoundingPMF);

                    return;
                }
            }
        }

        // Get the new success rate
        SuccessRate = (double)(Array - OldArray) / (double)SampleCount;

        if (SuccessRate < _RNG_MONTECARLO_MINSUCCESS)
            SuccessRate = _RNG_MONTECARLO_MINSUCCESS;
    }
}

uint32_t RNG_Poisson(RNG_Seed Seed, double Lambda)
{
    return RNG_MonteCarloUInt(Seed, &_RNG_PoissonPMF, &_RNG_PoissonBoundingPMF, &_RNG_PoissonBoundingSampler, &Lambda, 1);
}

uint32_t *RNG_PoissonArray(RNG_Seed Seed, double Lambda, size_t Size)
{
    return RNG_MonteCarloUIntArray(Seed, &_RNG_PoissonPMFArrayM, &_RNG_PoissonBoundingPMFArrayM, &_RNG_PoissonBoundingSamplerArrayM, &Lambda, 1, Size);
}

void RNG_PoissonArrayM(RNG_Seed Seed, double Lambda, uint32_t *Array, size_t Size)
{
    return RNG_MonteCarloUIntArrayM(Seed, &_RNG_PoissonPMFArrayM, &_RNG_PoissonBoundingPMFArrayM, &_RNG_PoissonBoundingSamplerArrayM, &Lambda, 1, Array, Size);
}

double RNG_PoissonPMF(uint32_t n, double Lambda)
{
    return exp(log(Lambda) * (double)n - Lambda - lgamma((double)n + 1));
}

double *RNG_PoissonPMFArray(uint32_t *n, double Lambda, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    RNG_PoissonPMFArrayM(n, Lambda, Array, Size);

    return Array;
}

double *RNG_PoissonPMFArrayM(uint32_t *n, double Lambda, double *Array, size_t Size)
{
    // Calculate constants
    double lLambda = log(Lambda);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        *List = exp(lLambda * (double)(*n) - Lambda - lgamma((double)(*n) + 1));
}

double _RNG_PoissonPMF(uint32_t n, void *Params)
{

}

double *_RNG_PoissonPMFArray(uint32_t *n, void *Params, size_t Size)
{

}

double *_RNG_PoissonPMFArrayM(uint32_t *n, void *Params, double *Array, size_t Size)
{

}

uint32_t _RNG_PoissonBoundingSampler(RNG_Seed Seed, void *Params)
{

}

uint32_t *_RNG_PoissonBoundingSamplerArray(RNG_Seed Seed, void *Params, size_t Size)
{

}

void _RNG_PoissonBoundingSamplerArrayM(RNG_Seed Seed, void *Params, uint32_t *Array, size_t Size)
{

}

double _RNG_PoissonBoundingPMF(uint32_t n, void *Params)
{

}

double *_RNG_PoissonBoundingPMFArray(uint32_t *n, void *Params, size_t Size)
{

}

double *_RNG_PoissonBoundingPMFArrayM(uint32_t *n, void *Params, double *Array, size_t Size)
{

}

// From libit
double _RNG_erfinv(double x)
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
        r = x * (((_RNG_ERFINV_A3 * x2 + _RNG_ERFINV_A2) * x2 + _RNG_ERFINV_A1) * x2 + _RNG_ERFINV_A0);
        r /= (((_RNG_ERFINV_B4 * x2 + _RNG_ERFINV_B3) * x2 + _RNG_ERFINV_B2) * x2 + _RNG_ERFINV_B1) * x2 + _RNG_ERFINV_B0;
    }
  
    else 
    {
        y = sqrt (-log ((1 - x) / 2));
        r = (((_RNG_ERFINV_C3 * y + _RNG_ERFINV_C2) * y + _RNG_ERFINV_C1) * y + _RNG_ERFINV_C0);
        r /= ((_RNG_ERFINV_D2 * y + _RNG_ERFINV_D1) * y + _RNG_ERFINV_D0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= _RNG_SQRTPI * (erf(r) - x) / (2 * exp(-r * r));
    r -= _RNG_SQRTPI * (erf(r) - x) / (2 * exp(-r * r));

    return r;
}

#endif // RANDOM2_H_INCLUDED
