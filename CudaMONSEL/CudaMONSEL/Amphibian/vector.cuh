#ifndef _VECTOR_CUH_
#define _VECTOR_CUH_

#include <cuda_runtime.h>

#include <stdio.h>

namespace amp
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static const int VECTOR_INITIAL_SIZE = 23;
#else
   static const int VECTOR_INITIAL_SIZE = 23;
#endif

   template<typename T>
   class vector
   {
   public:
      __host__ __device__ vector();
      __host__ __device__ ~vector();

   private:
      T** vec;
   };

   template<typename T>
   __host__ __device__ vector<T>::vector()
   {
      vec = new T*[VECTOR_INITIAL_SIZE];
      memset(vec, 0, sizeof(T*) * VECTOR_INITIAL_SIZE);
   }

   template<typename T>
   __host__ __device__ vector<T>::~vector()
   {

   }
}

#endif