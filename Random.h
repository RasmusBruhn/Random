#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define ERR_PREFIX RNG
#include <Error2.h>

typedef uint64_t RNG_Seed;

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

// Constants
#define RNG_MAX 18446744073709551616.
#define _RNG_MULTIPLIER 6364136223846793005
#define _RNG_CONSTANT 1442695040888963407

// Get a random uint64_t and update the seed to be that number
// Seed (RNG_Seed): The seed to use and update
#define RNG_IntFast(Seed) ((*Seed) = (*Seed) * _RNG_MULTIPLIER + _RNG_CONSTANT)

// Get a random double between 0 inclusive and 1 exclusive and update seed
// Seed (RNG_Seed): The seed to use and update
#define RNG_FloatFast(Seed) ((double) RNG_IntFast(Seed) / RNG_MAX)

#endif
