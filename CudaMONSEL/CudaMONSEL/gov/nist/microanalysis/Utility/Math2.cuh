#ifndef _MATH_2_CUH_
#define _MATH_2_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <cstdlib> 
#include <ctime>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

namespace Math2
{
   __host__ __device__ extern void minus3d(const double a[], const double b[], double res[]);
   __host__ __device__ extern void divide3d(const double a[], double b, double res[]);
   __host__ __device__ extern void plus3d(const double a[], const double b[], double res[]);
   __host__ __device__ extern void normalize3d(const double p[], double[]);
   __host__ __device__ extern void multiply3d(double a, const double b[], double[]);
   __host__ __device__ extern double distance3d(const double p1[], const double p2[]);
   __host__ __device__ extern double dot3d(const double a[], const double b[]);
   __host__ __device__ extern void cross3d(const double* a, const double* b, double res[]);
   __host__ __device__ extern double sqr(double x);
   __host__ __device__ extern double magnitude3d(const double p[]);
   __host__ __device__ extern void pointBetween3d(const double a[], const double b[], double f, double res[]);
   __host__ __device__ extern double toRadians(double deg);

   __host__ extern double random();
   __device__ extern double random(curandState&);
   __host__ extern int randomInt(int mod);
   __device__ extern int randomInt(int mod, curandState&);
   __host__ extern double expRand();
   __device__ extern double expRand(curandState&);
   __host__ extern double generateGaussianNoise(const double mean, const double stdDev);
   __device__ extern double generateGaussianNoise(const double mean, const double stdDev, curandState&);

   extern const double PI;
   extern const double ORIGIN_3D[];
   extern const double ONE[];
   extern const double X_AXIS[];
   extern const double Y_AXIS[];
   extern const double Z_AXIS[];
   extern const double MINUS_X_AXIS[];
   extern const double MINUS_Y_AXIS[];
   extern const double MINUS_Z_AXIS[];
}

#endif