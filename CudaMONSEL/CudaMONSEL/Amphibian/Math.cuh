#ifndef _MATH_CUH_
#define _MATH_CUH_

#include "cuda_runtime.h"

namespace Math
{
   __device__ __host__ double sqrt(double, double);
   __device__ __host__ double abs(double);

   __device__ __host__ double signum(double);
};
#endif
