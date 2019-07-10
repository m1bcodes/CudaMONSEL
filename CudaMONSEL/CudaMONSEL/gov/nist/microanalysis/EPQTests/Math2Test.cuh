#ifndef _MATH_2_TEST_CUH_
#define _MATH_2_TEST_CUH_

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

namespace Math2Test
{
   extern void testRandom1();
   extern void testRandom2();

   __device__ extern void testRandom1Cuda();
   __device__ extern void testRandom2Cuda();
}

#endif