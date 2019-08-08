#include "Amphibian\random.cuh"

#include <curand.h>
#include <curand_kernel.h>
#include <thread>
#include <random>
#include <time.h>

#if defined (_MSC_VER)  // Visual studio
#define thread_local __declspec(thread)
#elif defined (__GCC__) // GCC
#define thread_local __thread
#endif

namespace Random
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const float PI = 3.14159265358979323846f;
#else
   const float PI = 3.14159265358979323846f;
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

   // https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
   //Thread-safe function that returns a random number between min and max (inclusive).
   //This function takes ~142% the time that calling rand() would take. For this extra
   //cost you get a better uniform distribution and thread-safety.
   static int intRand(const int & min, const int & max)
   {
      static thread_local std::mt19937* generator = nullptr;
      if (!generator) generator = new std::mt19937(clock() + std::this_thread::get_id().hash());
      std::uniform_int_distribution<int> distribution(min, max);
      return distribution(*generator);
   }

   __host__ __device__ float random()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      // curand_uniform() can return between 0.0 (exclusive) and 1.0 (inclusive).
      int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      return 1 - curand_uniform(&states[threadId]); // excludes 1.0 but includes 0.0.
#else
      //return (double)rand() / RAND_MAX;
      // http://c-faq.com/lib/randrange.html
      // like Java, does not include 1
      //return (float)rand() / (RAND_MAX + 1);
      return (float)intRand(0, RAND_MAX) / (RAND_MAX + 1);
#endif
   }

   __host__ __device__ int randomInt(int mod)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      return (int)truncf((1 - curand_uniform(&states[threadId])) * mod);
#else
      //return rand() % mod;
      return intRand(0, mod-1);
#endif
   }

   __host__ __device__ float expRand()
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      float r = curand_uniform(&states[threadId]); // excludes 0.0 but includes 1.0.
      if (r >= 1) r = 1 - r;
      return -::log(r);
#else
      //float r = (float)rand() / RAND_MAX;
      //while (r <= 0 || r >= 1) r = (float)rand() / RAND_MAX;
      float r = intRand(1, RAND_MAX - 1);
      return -::logf(r/(float)RAND_MAX);
#endif
   }

   __host__ __device__ float generateGaussianNoise(const float mean, const float stdDev)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      int blockId = blockIdx.x + blockIdx.y * gridDim.x;
      int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
      return curand_normal(&states[threadId]) * stdDev + mean;
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
         u = (intRand(0, RAND_MAX) / ((float)RAND_MAX)) * 2.0f - 1.0f;
         v = (intRand(0, RAND_MAX) / ((float)RAND_MAX)) * 2.0f - 1.0f;
         s = u * u + v * v;
      } while ((s >= 1.0f) || (s == 0.0f));
      s = ::sqrtf(-2.0f * ::logf(s) / s);
      spare = v * s;
      return mean + stdDev * u * s;
#endif
   }
}