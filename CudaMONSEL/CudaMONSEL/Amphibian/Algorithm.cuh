#ifndef _ALGORITHM_CUH_
#define _ALGORITHM_CUH_

#include <cuda_runtime.h>

namespace Algorithm
{
   __host__ __device__ int binarySearch(const double a[], int l, int r, const double t);
   __host__ __device__ void quicksort(double a[], int l, int r);

   __host__ __device__ int binarySearch(const float a[], int l, int r, const float t);
}

#endif