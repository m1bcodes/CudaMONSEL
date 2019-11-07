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
   extern __host__ __device__ double fabs(double x);
   extern __host__ __device__ void minus3d(const double a[], const double b[], double res[]);
   extern __host__ __device__ void divide3d(const double a[], double b, double res[]);
   extern __host__ __device__ void plus3d(const double a[], const double b[], double res[]);
   extern __host__ __device__ void normalize3d(const double p[], double[]);
   extern __host__ __device__ void multiply3d(double a, const double b[], double[]);
   extern __host__ __device__ double distance3d(const double p1[], const double p2[]);
   extern __host__ __device__ double dot3d(const double a[], const double b[]);
   extern __host__ __device__ void cross3d(const double* a, const double* b, double res[]);
   extern __host__ __device__ double sqr(double x);
   extern __host__ __device__ double magnitude3d(const double p[]);
   extern __host__ __device__ void pointBetween3d(const double a[], const double b[], double f, double res[]);
   extern __host__ __device__ double toRadians(const double deg);

   extern __host__ __device__ void plus3d(const float a[], const float b[], float res[]);
   extern __host__ __device__ float distance3d(const float p1[], const float p2[]);
   extern __host__ __device__ float sqr(float x);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   extern __constant__ const double PI;
#else
   extern const double PI;
#endif
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