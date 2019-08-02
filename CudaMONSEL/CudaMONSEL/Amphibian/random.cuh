#ifndef _RANDOM_CUH_
#define _RANDOM_CUH_

#include <cuda_runtime.h>

namespace Random
{
   extern __host__ __device__ float random();
   extern __host__ __device__ int randomInt(int mod);
   extern __host__ __device__ float expRand();
   extern __host__ __device__ float generateGaussianNoise(const float mean, const float stdDev);

   extern __global__ void initCudaStates(const unsigned int n);
   extern __global__ void destroyCudaState();

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   extern __constant__ const float PI;
#else
   extern const float PI;
#endif
}

#endif