#include "Amphibian\random.cuh"

#include <curand.h>
#include <curand_kernel.h>

namespace Random
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const double PI = 3.14159265358979323846;
#else
   const double PI = 3.14159265358979323846;
#endif

   __device__ curandState *states;

   __global__ void initCudaStates(const unsigned int n)
   {
      states = new curandState[n];
      for (int i = 0; i < n; ++i) {
         curand_init(0, i, i, &states[i]);
      }
   }

   __global__ void destroyCudaState()
   {
      delete[] states;
   }

   __host__ __device__ float random()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      // curand_uniform() can return between 0.0 (exclusive) and 1.0 (inclusive).
      return 1 - curand_uniform(&states[threadIdx.x + blockDim.x * blockIdx.x]); // excludes 1.0 but includes 0.0.
#else
      //return (double)rand() / RAND_MAX;
      // http://c-faq.com/lib/randrange.html
      // like Java, does not include 1
      return (float)rand() / (RAND_MAX + 1);
#endif
   }

   __host__ __device__ int randomInt(int mod)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return (int)truncf((1 - curand_uniform(&states[threadIdx.x + blockDim.x * blockIdx.x])) * mod);
#else
      return rand() % mod;
#endif
   }

   __host__ __device__ float expRand()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      float r = curand_uniform(&states[threadIdx.x + blockDim.x * blockIdx.x]); // excludes 0.0 but includes 1.0.
      if (r >= 1) r = 1 - r;
      return -::log(r);
#else
      float r = (float)rand() / RAND_MAX;
      while (r <= 0 || r >= 1) r = (float)rand() / RAND_MAX;
      return -::log(r);
#endif
   }

   __host__ __device__ float generateGaussianNoise(const float mean, const float stdDev)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return curand_normal(&states[threadIdx.x + blockDim.x * blockIdx.x]) * stdDev + mean;
#else
      static bool hasSpare = false;
      static float spare;

      if (hasSpare) {
         hasSpare = false;
         return mean + stdDev * spare;
      }

      hasSpare = true;
      static float u, v, s;
      do {
         u = (rand() / ((float)RAND_MAX)) * 2.0 - 1.0;
         v = (rand() / ((float)RAND_MAX)) * 2.0 - 1.0;
         s = u * u + v * v;
      } while ((s >= 1.0) || (s == 0.0));
      s = ::sqrt(-2.0 * ::log(s) / s);
      spare = v * s;
      return mean + stdDev * u * s;
#endif
   }
}