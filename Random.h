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

typedef uint64_t RNG_Seed;
typedef struct ___RNG_BinomialParams _RNG_BinomialParams;

struct ___RNG_BinomialParams
{
    uint64_t N;
    double p;
};

// Generates a seed, the sed must be destroyed when no longer used
RNG_Seed *RNG_SeedGenerate();

// Creates a seed from a number, the sed must be destroyed when no longer used
// Value: The number to convert to a seed
RNG_Seed *RNG_SeedCreate(uint64_t Value);

// Destroys a seed, it must not be used after running this
void RNG_SeedDestroy(RNG_Seed *Seed);

// Get a uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
uint64_t RNG_Int(RNG_Seed *Seed, uint64_t Min, uint64_t Max);

// Get an array of uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
uint64_t *RNG_IntArray(RNG_Seed *Seed, uint64_t Min, uint64_t Max, size_t Size);

// Fill an array with uniform random int
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
uint64_t *RNG_IntArrayM(RNG_Seed *Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size);

// Get a uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
double RNG_Float(RNG_Seed *Seed, double Min, double Max);

// Get an array of uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Size: The size of the array
double *RNG_FloatArray(RNG_Seed *Seed, double Min, double Max, size_t Size);

// Fill an array with uniform random float
// Seed: The seed to use and update, NULL to use global seed
// Min: The minimum value
// Max: The maximum value
// Array: The array to fill
// Size: The size of the array
double *RNG_FloatArrayM(RNG_Seed *Seed, double Min, double Max, double *Array, size_t Size);

// Get a random number from an exponential distibution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
double RNG_Exp(RNG_Seed *Seed, double x0, double l);

// Get an array of random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
// Size: The size of the array
double *RNG_ExpArray(RNG_Seed *Seed, double x0, double l, size_t Size);

// Fill an array with random numbers from an exponential distribution (1 / l * exp(-x / l))
// Seed: The seed to use and update, NULL to use global seed
// x0: The minimu x possible
// l: The length scale of the distribution
// Array: The array to fill
// Size: The size of the array
double *RNG_ExpArrayM(RNG_Seed *Seed, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpPDFArrayM(double *x, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpCDFArrayM(double *x, double x0, double l, double *Array, size_t Size);

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
double *RNG_ExpICDFArrayM(double *y, double x0, double l, double *Array, size_t Size);

// Get a random number from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
double RNG_Normal(RNG_Seed *Seed, double Mu, double Sigma);

// Get an array of random numbers from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Size: The size of the array
double *RNG_NormalArray(RNG_Seed *Seed, double Mu, double Sigma, size_t Size);

// Fill an array with random numbers from a normal distribution (1 / (sigma * sqrt(2 * pi)) * exp(- (x - mu) ** 2 / (2 * sigma ** 2)))
// Seed: The seed to use and update, NULL to use global seed
// Mu: The mean of the distribution
// Sigma: The standard deviation of the distribution
// Array: The array to fill
// Size: The size of the array
double *RNG_NormalArrayM(RNG_Seed *Seed, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalPDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalCDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size);

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
double *RNG_NormalICDFArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size);

// Get a random number from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
uint64_t RNG_Poisson(RNG_Seed *Seed, double Lambda);

// Get an array of random numbers from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
// Size: The size of the array
uint64_t *RNG_PoissonArray(RNG_Seed *Seed, double Lambda, size_t Size);

// Fills an array of random numbers from the poisson distribution (lambda ** n * exp(lambda) / n!)
// Seed: The seed to use and update, NULL to use global seed
// Lambda: The mean value of the distribution
// Array: The array to fill
// Size: The size of the array
uint64_t *RNG_PoissonArrayM(RNG_Seed *Seed, double Lambda, uint64_t *Array, size_t Size);

// Get the PMF for a poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The position to get the PMF for
// Lambda: The mean value of the distribution
double RNG_PoissonPMF(uint64_t n, double Lambda);

// Get an array of PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Lambda: The mean value of the distribution
// Size: The size of the array
double *RNG_PoissonPMFArray(uint64_t *n, double Lambda, size_t Size);

// Fills an array of PMF values for the poisson distribution (lambda ** n * exp(lambda) / n!)
// n: The positions to get the PMF for
// Lambda: The mean value of the distribution
// Array: The array to fill
// Size: The size of the array
double *RNG_PoissonPMFArrayM(uint64_t *n, double Lambda, double *Array, size_t Size);

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

// Get a random number from the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// Seed: The seed to use and update, NULL to use global seed
// N: The number of events to draw
// p: The probability that an event is a success
uint64_t RNG_Binomial(RNG_Seed *Seed, uint64_t N, double p);

// Get an array of random numbers from the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// Seed: The seed to use and update, NULL to use global seed
// N: The number of events to draw
// p: The probability that an event is a success
// Size: The size of the array
uint64_t *RNG_BinomialArray(RNG_Seed *Seed, uint64_t N, double p, size_t Size);

// Fills an array of random numbers from the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// Seed: The seed to use and update, NULL to use global seed
// N: The number of events to draw
// p: The probability that an event is a success
// Array: The array to fill
// Size: The size of the array
uint64_t *RNG_BinomialArrayM(RNG_Seed *Seed, uint64_t N, double p, uint64_t *Array, size_t Size);

// Get the PMF for a binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// n: The position to get the PMF for
// N: The number of events to draw
// p: The probability that an event is a success
double RNG_BinomialPMF(uint64_t n, uint64_t N, double p);

// Get an array of PMF values for the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// n: The positions to get the PMF for
// N: The number of events to draw
// p: The probability that an event is a success
// Size: The size of the array
double *RNG_BinomialPMFArray(uint64_t *n, uint64_t N, double p, size_t Size);

// Fills an array of PMF values for the binomial distribution (N! / (n! * (N - n)!) p ** n * (1 - p) ** (N - n))
// n: The positions to get the PMF for
// N: The number of events to draw
// p: The probability that an event is a success
// Array: The array to fill
// Size: The size of the array
double *RNG_BinomialPMFArrayM(uint64_t *n, uint64_t N, double p, double *Array, size_t Size);

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

// Samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingPDFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
double RNG_MonteCarlo(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier);

// Creates an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
double *RNG_MonteCarloArray(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier, size_t Size);

// Fills an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PDFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPDFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PDF(x) <= BoundingMultiplier * BoundingPDF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
// Array: The array to fill
double *RNG_MonteCarloArrayM(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier, double *Array, size_t Size);

// Samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingPMFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position to get the PDF, Params: A pointer to some type containing all parameters for the PDF
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
uint64_t RNG_MonteCarloUInt(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier);

// Creates an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPMFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
uint64_t *RNG_MonteCarloUIntArray(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier, size_t Size);

// Fills an array with samples from some distribution using Monte Carlo simulation
// Seed: The seed to use and update, NULL to use global seed
// PMFArrayM: The PDF function for the distribution to sample from, it does not need to be normalized, takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingPMFArrayM: The PDF for the bounding function, it does not need to normalized but for all x: PMF(x) <= BoundingMultiplier * BoundingPMF(x), takes the arguments: x: The position array to get the PDF, Params: A pointer to some type containing all parameters for the PDF, Array: The array to fill, Size: The size of the array
// BoundingSamplerArrayM: The function to get a sample form the bounding function, takes the arguments: Seed: The seed to use, Params: A pointer to some type containing all parameters for the sampler, Array: The array to fill, Size: The size of the array
// Params: A pointer to some type containing all parameters for the PDFs
// BoundingMultiplier: The number to multiply the BoundingPDF with to allow it to be larger than the PDF
// Size: The size of the array
// Array: The array to fill
uint64_t *RNG_MonteCarloUIntArrayM(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *x, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *x, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier, uint64_t *Array, size_t Size);

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
#define _RNG_1_E 0.36787944117144233
#define _RNG_E 2.718281828459045

RNG_Seed _RNG_GlobalSeed = 0;

// Get a random uint64_t and update the seed to be that number
// Seed (RNG_Seed): The seed to use and update
#define RNG_IntFast(Seed) ((*Seed) = (*Seed) * _RNG_MULTIPLIER + _RNG_CONSTANT)

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (RNG_Seed): The seed to use and update
#define RNG_FloatFast(Seed) ((double) RNG_IntFast(Seed) / RNG_MAX)

RNG_Seed *RNG_SeedGenerate()
{
    // Get the seed
    uint64_t Value = (uint64_t) time(NULL) ^ (uint64_t) clock();

    // Create seed
    RNG_Seed *Seed = RNG_SeedCreate(Value);

    if (Seed == NULL)
        _RNG_ErrorAdd(_RNG_ERRORMES_CREATESEED);

    return Seed;
}

RNG_Seed *RNG_SeedCreate(uint64_t Value)
{
    // Allocate memory for it
    RNG_Seed *Seed = (RNG_Seed *)malloc(sizeof(RNG_Seed));

    if (Seed == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t));
        return NULL;
    }

    // Set the value
    *Seed = Value;

    return Seed;
}

void RNG_SeedDestroy(RNG_Seed *Seed)
{
    free(Seed);
}

uint64_t RNG_Int(RNG_Seed *Seed, uint64_t Min, uint64_t Max)
{
    uint64_t Value;

    if (RNG_IntArrayM(Seed, Min, Max, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "int");
        return -1;
    }

    return Value;
}

uint64_t *RNG_IntArray(RNG_Seed *Seed, uint64_t Min, uint64_t Max, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    if (RNG_IntArrayM(Seed, Min, Max, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "int");
        free(Array);
        return NULL;
    }

    return Array;
}

uint64_t *RNG_IntArrayM(RNG_Seed *Seed, uint64_t Min, uint64_t Max, uint64_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Make sure max is not small than min
    if (Max < Min)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_LARGERTHAN, "Max", "Min", Max, Min);
        return NULL;
    }

    // Fill memory
    for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + (RNG_IntFast(Seed) % (1 + Max - Min));

    return Array;
}

double RNG_Float(RNG_Seed *Seed, double Min, double Max)
{
    double Value;

    if (RNG_FloatArrayM(Seed, Min, Max, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "float");
        return NAN;
    }

    return Value;
}

double *RNG_FloatArray(RNG_Seed *Seed, double Min, double Max, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_FloatArrayM(Seed, Min, Max, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "float");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_FloatArrayM(RNG_Seed *Seed, double Min, double Max, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    if (Max < Min)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_LARGERTHANF, "Max", "Min", Max, Min);
        return NULL;
    }

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        *List = Min + RNG_FloatFast(Seed) * (Max - Min);

    return Array;
}

double RNG_Exp(RNG_Seed *Seed, double x0, double l)
{
    double Value;

    if (RNG_ExpArrayM(Seed, x0, l, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "exp");
        return NAN;
    }

    return Value;
}

double *RNG_ExpArray(RNG_Seed *Seed, double x0, double l, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_ExpArrayM(Seed, x0, l, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "exp");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_ExpArrayM(RNG_Seed *Seed, double x0, double l, double *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    if (l < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "l", l);
        return NULL;
    }

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FloatFast(Seed);

        // Get from distribution
        *List = x0 - l * log(1 - Uniform);
    }

    return Array;
}

double RNG_ExpPDF(double x, double x0, double l)
{
    double Value;

    if (RNG_ExpPDFArrayM(&x, x0, l, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PDF, "exp");
        return NAN;
    }

    return Value;
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
    if (RNG_ExpPDFArrayM(x, x0, l, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PDF, "exp");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_ExpPDFArrayM(double *x, double x0, double l, double *Array, size_t Size)
{
    if (l < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "l", l);
        return NULL;
    }

    if (l == 0)
    {
        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
            *List = ((*x == x0) ? (INFINITY) : (0));

        return Array;
    }

    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < x0) ? (0) : (c * exp(-c * (*x - x0))));

    return Array;
}

double RNG_ExpCDF(double x, double x0, double l)
{
    double Value;

    if (RNG_ExpCDFArrayM(&x, x0, l, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_CDF, "exp");
        return NAN;
    }

    return Value;
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
    if (RNG_ExpCDFArrayM(x, x0, l, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_CDF, "exp");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_ExpCDFArrayM(double *x, double x0, double l, double *Array, size_t Size)
{
    if (l < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "l", l);
        return NULL;
    }

    if (l == 0)
    {
        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
            *List = ((*x < x0) ? (0) : (1));

        return Array;
    }

    // Calculate constants
    double c = 1 / l;

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = ((*x < x0) ? (0) : (1 - exp(-c * (*x - x0))));

    return Array;
}

double RNG_ExpICDF(double y, double x0, double l)
{
    double Value;

    if (RNG_ExpICDFArrayM(&y, x0, l, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_ICDF, "exp");
        return NAN;
    }

    return Value;
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
    if (RNG_ExpICDFArrayM(y, x0, l, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_ICDF, "exp");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_ExpICDFArrayM(double *y, double x0, double l, double *Array, size_t Size)
{
    if (l < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "l", l);
        return NULL;
    }

    if (l == 0)
    {
        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
            *List = x0;

        return Array;
    }

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = ((*y < 0) ? (x0) : ((*y >= 1) ? (INFINITY) : (x0 - l * log(1 - *y))));

    return Array;
}

double RNG_Normal(RNG_Seed *Seed, double Mu, double Sigma)
{
    double Value;

    if (RNG_NormalArrayM(Seed, Mu, Sigma, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "normal");
        return NAN;
    }

    return Value;
}

double *RNG_NormalArray(RNG_Seed *Seed, double Mu, double Sigma, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_NormalArrayM(Seed, Mu, Sigma, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "normal");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_NormalArrayM(RNG_Seed *Seed, double Mu, double Sigma, double *Array, size_t Size)
{
    if (Sigma < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Sigma", Sigma);
        return NULL;
    }

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

    return Array;
}

double RNG_NormalPDF(double x, double Mu, double Sigma)
{
    double Value;

    if (RNG_NormalPDFArrayM(&x, Mu, Sigma, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PDF, "normal");
        return NAN;
    }

    return Value;
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
    if (RNG_NormalPDFArrayM(x, Mu, Sigma, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PDF, "normal");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_NormalPDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    if (Sigma < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Sigma", Sigma);
        return NULL;
    }

    if (Sigma == 0)
    {
        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
            *List = ((*x == Mu) ? (INFINITY) : (0));

        return Array;
    }

    // Calculate constants
    double A = 1 / (Sigma * _RNG_SQRTPI * _RNG_SQRT2);
    double B = 0.5 / (Sigma * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = A * exp(-B * (*x - Mu) * (*x - Mu));

    return Array;
}

double RNG_NormalCDF(double x, double Mu, double Sigma)
{
    double Value;

    if (RNG_NormalCDFArrayM(&x, Mu, Sigma, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_CDF, "normal");
        return NAN;
    }

    return Value;
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
    if (RNG_NormalCDFArrayM(x, Mu, Sigma, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_CDF, "normal");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_NormalCDFArrayM(double *x, double Mu, double Sigma, double *Array, size_t Size)
{
    if (Sigma < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Sigma", Sigma);
        return NULL;
    }

    if (Sigma == 0)
    {
        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
            *List = ((*x < Mu) ? (0) : (1));

        return Array;
    }

    // Calculate constants
    double A = 1 / (_RNG_SQRT2 * Sigma);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++x)
        *List = 0.5 * (1 + erf(A * (*x - Mu)));

    return Array;
}

double RNG_NormalICDF(double y, double Mu, double Sigma)
{
    double Value;

    if (RNG_NormalICDFArrayM(&y, Mu, Sigma, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_ICDF, "normal");
        return NAN;
    }

    return Value;
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
    if (RNG_NormalICDFArrayM(y, Mu, Sigma, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_ICDF, "normal");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_NormalICDFArrayM(double *y, double Mu, double Sigma, double *Array, size_t Size)
{
    if (Sigma < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Sigma", Sigma);
        return NULL;
    }

    // Calculate constants
    double A = 2 * Sigma;
 
    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++y)
        *List = ((*y < 0) ? (-INFINITY) : ((*y >= 1) ? (INFINITY) : (Mu + A * _RNG_erfinv(2 * *y - 1))));

    return Array;
}

uint64_t RNG_Poisson(RNG_Seed *Seed, double Lambda)
{
    uint64_t Value;

    if (RNG_PoissonArrayM(Seed, Lambda, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "poisson");
        return -1;
    }

    return Value;
}

uint64_t *RNG_PoissonArray(RNG_Seed *Seed, double Lambda, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_PoissonArrayM(Seed, Lambda, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "poisson");
        free(Array);
        return NULL;
    }

    return Array;
}

uint64_t *RNG_PoissonArrayM(RNG_Seed *Seed, double Lambda, uint64_t *Array, size_t Size)
{
    if (Lambda < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Lambda", Lambda);
        return NULL;
    }

    if (Lambda == 0)
    {
        for (uint64_t *List = Array, *EndList = Array + Size; List < EndList; ++List)
            *List = 0;

        return Array;
    }

    if (RNG_MonteCarloUIntArrayM(Seed, &_RNG_PoissonPMFArrayM, &_RNG_PoissonBoundingPMFArrayM, &_RNG_PoissonBoundingSamplerArrayM, &Lambda, 1, Array, Size) == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_GENERATE, "poisson");
        return NULL;
    }

    return Array;
}

double RNG_PoissonPMF(uint64_t n, double Lambda)
{
    double Value;

    if (RNG_PoissonPMFArrayM(&n, Lambda, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "poisson");
        return NAN;
    }

    return Value;
}

double *RNG_PoissonPMFArray(uint64_t *n, double Lambda, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_PoissonPMFArrayM(n, Lambda, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "poisson");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_PoissonPMFArrayM(uint64_t *n, double Lambda, double *Array, size_t Size)
{
    if (Lambda < 0)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_NOTNEGATIVE, "Lambda", Lambda);
        return NULL;
    }

    if (Lambda == 0)
    {
        for (double *List = Array, *EndList = Array + Size; List < EndList; ++List, ++n)
            *List = ((*n == 0) ? (1) : (0));

        return Array;
    }

    if (_RNG_PoissonPMFArrayM(n, &Lambda, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "poisson");
        return NULL;
    }

    return Array;
}

double *_RNG_PoissonPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size)
{
    double Lambda = *(double *)Params;

    // Calculate constants
    double lLambda = log(Lambda);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        *List = exp(lLambda * (double)(*n) - Lambda - lgamma((double)(*n) + 1));

    return Array;
}

uint64_t *_RNG_PoissonBoundingSamplerArrayM(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size)
{
    // Unpack params
    double Lambda = *(double *)Params;

    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // If lambda is too small
    if (Lambda < 1)
    {
        for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        {
            // Get uniform number
            double Uniform = RNG_FloatFast(Seed);

            // Get from distribution
            *List = (uint64_t)(- log(1 - Uniform));
        }

        return Array;
    }

    // Calulate constants
    double LogLambda = log(Lambda);
    double SqrtLambda = sqrt(Lambda);
    double nm = ceil(Lambda * exp(-1 / SqrtLambda)) - 1;
    double np = ceil(Lambda * exp(1 / SqrtLambda)) - 1;
    double A = exp(LogLambda * (nm - 0.5) - Lambda - (nm - Lambda) / SqrtLambda - lgamma(nm + 1)) / (exp(1 / SqrtLambda) - 1);
    double B = exp(LogLambda * (np - 0.5) - Lambda + (np - Lambda) / SqrtLambda - lgamma(np + 1)) / (1 - exp(-1 / SqrtLambda));
    double PA = A * (1 - exp(-SqrtLambda)) / (B + A * (1 - exp(-SqrtLambda)));
    double AmpA = (1 - exp(-SqrtLambda));

    // Fill memory
    for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
    {
        // Get uniform number
        double Uniform = RNG_FloatFast(Seed);

        // It is part of the first exp
        if (Uniform < PA)
        {
            Uniform /= PA;
            *List = (uint64_t)(Lambda + SqrtLambda * log(1 - AmpA * Uniform));
        }

        // It is part of the second exp
        else
        {
            Uniform = (Uniform - PA) / (1 - PA);
            *List = (uint64_t)(Lambda - SqrtLambda * log(1 - Uniform));
        }
    }

    return Array;
}

double *_RNG_PoissonBoundingPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size)
{
    // Unpack params
    double Lambda = *(double *)Params;

    // If lambda is too small
    if (Lambda < 1)
    {
        double Amp = 0.5 * _RNG_E;

        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
            *List = Amp * exp(-(double)(*n));

        return Array;
    }

    // Calulate constants
    double LogLambda = log(Lambda);
    double SqrtLambda = sqrt(Lambda);
    double Std = 1 / SqrtLambda;
    double nm = ceil(Lambda * exp(-1 / SqrtLambda)) - 1;
    double np = ceil(Lambda * exp(1 / SqrtLambda)) - 1;
    double AmpA = exp(LogLambda * nm - Lambda - (nm - Lambda) / SqrtLambda - lgamma(nm + 1));
    double AmpB = exp(LogLambda * np - Lambda + (np - Lambda) / SqrtLambda - lgamma(np + 1));

    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
    {
        uint64_t nThreshold = (uint64_t)Lambda;

        // The first exp
        if (*n < nThreshold)
            *List = AmpA * exp(Std * ((double)(*n) - Lambda));

        // The second exp
        else if (*n > nThreshold)
            *List = AmpB * exp(-Std * ((double)(*n) - Lambda));

        // At the threshold
        else
            *List = AmpA * (1 - exp(Std * ((double)nThreshold - Lambda))) / (exp(Std) - 1) + AmpB * (1 - exp(-Std * ((double)nThreshold + 1 - Lambda))) / (1 - exp(-Std));
    }

    return Array;
}

uint64_t RNG_Binomial(RNG_Seed *Seed, uint64_t N, double p)
{
    uint64_t Value;

    if (RNG_BinomialArrayM(Seed, N, p, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "binomial");
        return -1;
    }

    return Value;
}

uint64_t *RNG_BinomialArray(RNG_Seed *Seed, uint64_t N, double p, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_BinomialArrayM(Seed, N, p, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "binomial");
        free(Array);
        return NULL;
    }

    return Array;
}

uint64_t *RNG_BinomialArrayM(RNG_Seed *Seed, uint64_t N, double p, uint64_t *Array, size_t Size)
{
    if (p < 0 || p > 1)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_INRANGE, "p", 0, 1, p);
        return NULL;
    }

    if (p == 0 || N == 0)
    {
        for (uint64_t *List = Array, *EndList = Array + Size; List < EndList; ++List)
            *List = 0;

        return Array;
    }

    if (p == 1)
    {
        for (uint64_t *List = Array, *EndList = Array + Size; List < EndList; ++List)
            *List = N;

        return Array;
    }

    if (N == 1)
    {
        // Get the global seed
        extern uint64_t _RNG_GlobalSeed;

        if (Seed == NULL)
            Seed = &_RNG_GlobalSeed;

        for (uint64_t *List = Array, *EndList = Array + Size; List < EndList; ++List)
        {
            double Uniform = RNG_FloatFast(Seed);

            if (Uniform < p)
                *List = 0;

            else
                *List = 1;
        }

        return Array;
    }

    _RNG_BinomialParams Params = {.N = N, .p = p};

    if (RNG_MonteCarloUIntArrayM(Seed, &_RNG_BinomialPMFArrayM, &_RNG_BinomialBoundingPMFArrayM, &_RNG_BinomialBoundingSamplerArrayM, &Params, 1, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "binomial");
        return NULL;
    }

    return Array;
}

double RNG_BinomialPMF(uint64_t n, uint64_t N, double p)
{
    double Value;

    if (RNG_BinomialPMFArrayM(&n, N, p, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "binomial");
        return NAN;
    }

    return Value;
}

double *RNG_BinomialPMFArray(uint64_t *n, uint64_t N, double p, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_BinomialPMFArrayM(n, N, p, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "binomial");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_BinomialPMFArrayM(uint64_t *n, uint64_t N, double p, double *Array, size_t Size)
{
    if (p < 0 || p > 1)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_INRANGE, "p", 0, 1, p);
        return NULL;
    }

    if (p == 0)
    {
        for (double *List = Array, *EndList = Array + Size; List < EndList; ++List, ++n)
            *List = ((*n == 0) ? (1) : (0));

        return Array;
    }

    if (p == 1)
    {
        for (double *List = Array, *EndList = Array + Size; List < EndList; ++List, ++n)
            *List = ((*n == N) ? (1) : (0));

        return Array;
    }

    _RNG_BinomialParams Params = {.N = N, .p = p};

    if (_RNG_BinomialPMFArrayM(n, &Params, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_PMF, "binomial");
        return NULL;
    }
    
    return Array;
}

double *_RNG_BinomialPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size)
{
    double N = (double)((_RNG_BinomialParams *)Params)->N;
    double p = ((_RNG_BinomialParams *)Params)->p;

    // Calculate constants
    double lgammaN = lgamma(N + 1);
    double logP = log(p);
    double logQ = log(1 - p);

    // Fill memory
    for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        *List = ((*n > N) ? (0) : (exp(lgammaN - lgamma((double)(*n) + 1) - lgamma(N - (double)(*n) + 1) + (double)(*n) * logP + (N - (double)(*n)) * logQ)));

    return Array;
}

uint64_t *_RNG_BinomialBoundingSamplerArrayM(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Unpack params
    double N = (double)((_RNG_BinomialParams *)Params)->N;
    double p = ((_RNG_BinomialParams *)Params)->p;

    // For low p
    if (N * p < 0.5)
    {
        double DecayRate = 1 / sqrt((1 - 1 / (2 * N)) / 2);
        double DecayLength = 1 / DecayRate;
        double mu = 1;
        double np = ceil((N * p - exp(-DecayRate) * (1 - p)) / (p + exp(-DecayRate) * (1 - p)));
        double IntB = exp(lgamma(N + 1) - lgamma(np + 1) - lgamma(N - np + 1) + DecayRate * (np - mu) + log(p) * np + log(1 - p) * (N - np)) * (1 - exp(-DecayRate * (N + 1 - mu))) / (1 - exp(-DecayRate));
        double LIntB = 1 - exp(-DecayRate * (N + 1 - mu));
        double P0 = exp(log(1 - p) * N);
        double PB = IntB / (P0 + IntB);

        // Fill memory
        for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        {
            // Get uniform number
            double Uniform = RNG_FloatFast(Seed);

            // It is part of the exp
            if (Uniform < PB)
            {
                Uniform /= PB;
                *List = (uint64_t)(mu - DecayLength * log(1 - LIntB * Uniform));
            }

            // It is 0
            else
                *List = 0;
        }
    }

    // For high p
    else if (N * p > N - 0.5)
    {
        double DecayRate = 1 / sqrt((1 - 1 / (2 * N)) / 2);
        double DecayLength = 1 / DecayRate;
        double mu = N;
        double nm = ceil((N * p - exp(DecayRate) * (1 - p)) / (p + exp(DecayRate) * (1 - p)));
        double IntA = exp(lgamma(N + 1) - lgamma(nm + 1) - lgamma(N - nm + 1) - DecayRate * (nm - mu) + log(p) * nm + log(1 - p) * (N - nm)) * (1 - exp(-DecayRate * mu)) / (exp(DecayRate) - 1);
        double LIntA = 1 - exp(-DecayRate * mu);
        double P0 = exp(log(p) * N);
        double PA = IntA / (P0 + IntA);

        // Fill memory
        for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        {
            // Get uniform number
            double Uniform = RNG_FloatFast(Seed);

            // It is part of the exp
            if (Uniform < PA)
            {
                Uniform /= PA;
                *List = (uint64_t)(mu + DecayLength * log(1 - LIntA * Uniform));
            }

            // It is N
            else
                *List = N;
        }
    }

    // For normal values of p
    else
    {
        double DecayRate = 1 / sqrt(N * p * (1 - p));
        double DecayLength = 1 / DecayRate;
        double mu = N * p + 0.5;
        double nm = ceil((N * p - exp(DecayRate) * (1 - p)) / (p + exp(DecayRate) * (1 - p)));
        double np = ceil((N * p - exp(-DecayRate) * (1 - p)) / (p + exp(-DecayRate) * (1 - p)));
        double IntA = exp(lgamma(N + 1) - lgamma(nm + 1) - lgamma(N - nm + 1) - DecayRate * (nm - mu) + log(p) * nm + log(1 - p) * (N - nm)) * (1 - exp(-DecayRate * mu)) / (exp(DecayRate) - 1);
        double IntB = exp(lgamma(N + 1) - lgamma(np + 1) - lgamma(N - np + 1) + DecayRate * (np - mu) + log(p) * np + log(1 - p) * (N - np)) * (1 - exp(-DecayRate * (N + 1 - mu))) / (1 - exp(-DecayRate));
        double LIntA = 1 - exp(-DecayRate * mu);
        double LIntB = 1 - exp(-DecayRate * (N + 1 - mu));
        double PA = IntA / (IntA + IntB);

        // Fill memory
        for (uint64_t *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List)
        {
            // Get uniform number
            double Uniform = RNG_FloatFast(Seed);

            // It is part of the first exp
            if (Uniform < PA)
            {
                Uniform /= PA;
                *List = (uint64_t)(mu + DecayLength * log(1 - LIntA * Uniform));
            }

            // It is part of the second exp
            else
            {
                Uniform = (Uniform - PA) / (1 - PA);
                *List = (uint64_t)(mu - DecayLength * log(1 - LIntB * Uniform));
            }
        }
    }

    return Array;
}

double *_RNG_BinomialBoundingPMFArrayM(uint64_t *n, const void *Params, double *Array, size_t Size)
{
    // Unpack params
    double N = (double)((_RNG_BinomialParams *)Params)->N;
    double p = ((_RNG_BinomialParams *)Params)->p;

    // It is small p
    if (N * p < 0.5)
    {
        double DecayRate = 1 / sqrt((1 - 1 / (2 * N)) / 2);
        double mu = 1;
        double np = ceil((N * p - exp(-DecayRate) * (1 - p)) / (p + exp(-DecayRate) * (1 - p)));
        double AmpB = exp(lgamma(N + 1) - lgamma(np + 1) - lgamma(N - np + 1) + DecayRate * (np - mu) + log(p) * np + log(1 - p) * (N - np));
        double P0 = exp(log(1 - p) * N);

        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        {
            if (*n == 0)
                *List = P0;

            else
                *List = AmpB * exp(-DecayRate * ((double)(*n) - mu));
        }
    }

    // For large p
    else if (N * p > N - 0.5)
    {
        double DecayRate = 1 / sqrt((1 - 1 / (2 * N)) / 2);
        double mu = N;
        double nm = ceil((N * p - exp(DecayRate) * (1 - p)) / (p + exp(DecayRate) * (1 - p)));
        double AmpA = exp(lgamma(N + 1) - lgamma(nm + 1) - lgamma(N - nm + 1) - DecayRate * (nm - mu) + log(p) * nm + log(1 - p) * (N - nm));
        double P0 = exp(log(p) * N);

        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        {
            if (*n == N)
                *List = P0;

            else
                *List = AmpA * exp(DecayRate * ((double)(*n) - mu));
        }
    }

    // Normal values for p
    else
    {
        double DecayRate = 1 / sqrt(N * p * (1 - p));
        double mu = N * p + 0.5;
        double nm = ceil((N * p - exp(DecayRate) * (1 - p)) / (p + exp(DecayRate) * (1 - p)));
        double np = ceil((N * p - exp(-DecayRate) * (1 - p)) / (p + exp(-DecayRate) * (1 - p)));
        double AmpA = exp(lgamma(N + 1) - lgamma(nm + 1) - lgamma(N - nm + 1) - DecayRate * (nm - mu) + log(p) * nm + log(1 - p) * (N - nm));
        double AmpB = exp(lgamma(N + 1) - lgamma(np + 1) - lgamma(N - np + 1) + DecayRate * (np - mu) + log(p) * np + log(1 - p) * (N - np));
        double ModA = AmpA / (exp(DecayRate) - 1);
        double ModB = AmpB / (1 - exp(-DecayRate));
        uint64_t muFloor = (uint64_t)floor(mu);
        double muLeft = mu - (double)muFloor;

        for (double *List = Array, *ListEnd = Array + Size; List < ListEnd; ++List, ++n)
        {
            if (*n < muFloor)
                *List = AmpA * exp(DecayRate * ((double)(*n) - mu));

            else if (*n > muFloor)
                *List = AmpB * exp(-DecayRate * ((double)(*n) - mu));

            else
                *List = ModA * (1 - exp(-DecayRate * muLeft)) + ModB * (1 - exp(-DecayRate * (1 - muLeft)));
        }
    }

    return Array;
}

double RNG_MonteCarlo(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier)
{
    double Value;

    if (RNG_MonteCarloArrayM(Seed, PDFArrayM, BoundingPDFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "Monte Carlo");
        return NAN;
    }

    return Value;
}

double *RNG_MonteCarloArray(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier, size_t Size)
{
    // Get memory
    double *Array = (double *)malloc(sizeof(double) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_MonteCarloArrayM(Seed, PDFArrayM, BoundingPDFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "Monte Carlo");
        free(Array);
        return NULL;
    }

    return Array;
}

double *RNG_MonteCarloArrayM(RNG_Seed *Seed, double *(*PDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingPDFArrayM)(double *x, const void *Params, double *Array, size_t Size), double *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, double *Array, size_t Size), const void *Params, double BoundingMultiplier, double *Array, size_t Size)
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
        return NULL;
    }

    double *PDF = (double *)malloc(sizeof(double) * Size);

    if (PDF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(x);
        return NULL;
    }

    double *BoundingPDF = (double *)malloc(sizeof(double) * Size);

    if (BoundingPDF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(x);
        free(PDF);
        return NULL;
    }

    // Set end point
    double *EndList = Array + Size;
    double *List = Array;
    double SuccessRate = 1;

    // Loop until done
    while (true)
    {
        // Get x samples
        size_t SampleCount = (size_t)((double)(EndList - List) / SuccessRate * _RNG_MONTECARLO_BUFFER);

        if (SampleCount > Size)
            SampleCount = Size;

        if (BoundingSamplerArrayM(Seed, Params, x, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOSAMPLE);
            free(x);
            free(PDF);
            free(BoundingPDF);
            return NULL;
        }

        // Get the PDF values
        if (PDFArrayM(x, Params, PDF, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOPDF);
            free(x);
            free(PDF);
            free(BoundingPDF);
            return NULL;
        }

        if (BoundingPDFArrayM(x, Params, BoundingPDF, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOBOUNDINGPDF);
            free(x);
            free(PDF);
            free(BoundingPDF);
            return NULL;
        }

        // Check if they are good
        double *OldList = List;

        for (double *xSample = x, *PDFSample = PDF, *BoundingPDFSample = BoundingPDF, *xSampleEnd = xSample + SampleCount; xSample < xSampleEnd; ++xSample, ++PDFSample, ++BoundingPDFSample)
        {
            // Get probability for success
            double Prob = *PDFSample / (BoundingMultiplier * *BoundingPDFSample);
            double Uniform = RNG_FloatFast(Seed);

            // Check if it should be kept
            if (Uniform < Prob)
            {
                *(List++) = *xSample;

                // Check if it is done
                if (List >= EndList)
                {
                    free(x);
                    free(PDF);
                    free(BoundingPDF);

                    return Array;
                }
            }
        }

        // Get the new success rate
        SuccessRate = (double)(List - OldList) / (double)SampleCount;

        if (SuccessRate < _RNG_MONTECARLO_MINSUCCESS)
            SuccessRate = _RNG_MONTECARLO_MINSUCCESS;
    }
}

uint64_t RNG_MonteCarloUInt(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier)
{
    uint64_t Value;

    if (RNG_MonteCarloUIntArrayM(Seed, PMFArrayM, BoundingPMFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, &Value, 1) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "Monte Carlo");
        return -1;
    }

    return Value;
}

uint64_t *RNG_MonteCarloUIntArray(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier, size_t Size)
{
    // Get memory
    uint64_t *Array = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (Array == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    // Get numbers
    if (RNG_MonteCarloUIntArrayM(Seed, PMFArrayM, BoundingPMFArrayM, BoundingSamplerArrayM, Params, BoundingMultiplier, Array, Size) == NULL)
    {
        _RNG_ErrorAdd(_RNG_ERRORMES_GENERATE, "Monte Carlo");
        free(Array);
        return NULL;
    }

    return Array;
}

uint64_t *RNG_MonteCarloUIntArrayM(RNG_Seed *Seed, double *(*PMFArrayM)(uint64_t *n, const void *Params, double *Array, size_t Size), double *(*BoundingPMFArrayM)(uint64_t *x, const void *Params, double *Array, size_t Size), uint64_t *(*BoundingSamplerArrayM)(RNG_Seed *Seed, const void *Params, uint64_t *Array, size_t Size), const void *Params, double BoundingMultiplier, uint64_t *Array, size_t Size)
{
    // Get the global seed
    extern uint64_t _RNG_GlobalSeed;

    if (Seed == NULL)
        Seed = &_RNG_GlobalSeed;

    // Get memory
    uint64_t *n = (uint64_t *)malloc(sizeof(uint64_t) * Size);

    if (n == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(uint64_t) * Size);
        return NULL;
    }

    double *PMF = (double *)malloc(sizeof(double) * Size);

    if (PMF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(n);
        return NULL;
    }

    double *BoundingPMF = (double *)malloc(sizeof(double) * Size);

    if (BoundingPMF == NULL)
    {
        _RNG_ErrorSet(_RNG_ERRORMES_MALLOC, sizeof(double) * Size);
        free(n);
        free(PMF);
        return NULL;
    }

    // Set end point
    uint64_t *EndList = Array + Size;
    uint64_t *List = Array;
    double SuccessRate = 1;

    // Loop until done
    while (true)
    {
        // Get x samples
        size_t SampleCount = (size_t)((double)(EndList - List) / SuccessRate * _RNG_MONTECARLO_BUFFER);

        if (SampleCount > Size)
            SampleCount = Size;

        if (BoundingSamplerArrayM(Seed, Params, n, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOSAMPLE);
            free(n);
            free(PMF);
            free(BoundingPMF);
            return NULL;
        }

        // Get the PDF values
        if (PMFArrayM(n, Params, PMF, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOPMF);
            free(n);
            free(PMF);
            free(BoundingPMF);
            return NULL;
        }

        if (BoundingPMFArrayM(n, Params, BoundingPMF, SampleCount) == NULL)
        {
            _RNG_ErrorAdd(_RNG_ERRORMES_MONTECARLOBOUNDINGPMF);
            free(n);
            free(PMF);
            free(BoundingPMF);
            return NULL;
        }

        // Check if they are good
        uint64_t *OldList = List;
        double *PMFSample = PMF;
        double *BoundingPMFSample = BoundingPMF;

        for (uint64_t *nSample = n, *nSampleEnd = nSample + SampleCount; nSample < nSampleEnd; ++nSample, ++PMFSample, ++BoundingPMFSample)
        {
            // Get probability for success
            double Prob = *PMFSample / (BoundingMultiplier * *BoundingPMFSample);
            double Uniform = RNG_FloatFast(Seed);

            // Check if it should be kept
            if (Uniform < Prob)
            {
                *(List++) = *nSample;

                // Check if it is done
                if (List >= EndList)
                {
                    free(n);
                    free(PMF);
                    free(BoundingPMF);

                    return Array;
                }
            }
        }

        // Get the new success rate
        SuccessRate = (double)(List - OldList) / (double)SampleCount;

        if (SuccessRate < _RNG_MONTECARLO_MINSUCCESS)
            SuccessRate = _RNG_MONTECARLO_MINSUCCESS;
    }
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
