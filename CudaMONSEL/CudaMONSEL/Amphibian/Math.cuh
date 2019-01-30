#ifndef _MATH_CUH_
#define _MATH_CUH_

#include "cuda_runtime.h"

class Math
{
public:
   __device__ __host__ static double sqrt(double, double);
   __device__ __host__ static double abs(double);
};
#endif
