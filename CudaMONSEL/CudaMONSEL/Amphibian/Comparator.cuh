#ifndef _COMPARATOR_CUH_
#define _COMPARATOR_CUH_

#include <cuda_runtime.h>

namespace Comparator
{
   typedef bool(*pIntComparator)(int, int);
   typedef bool(*pDoubleComparator)(double, double);
   typedef bool(*pLongComparator)(long, long);

   template<typename T>
   __host__ __device__ bool BuildCmp(T a, T b)
   {
      return a == b;
   }

   template<typename T>
   __host__ __device__ vector::vector()
   {

   }

   template<typename T>
   __host__ __device__ vector::~vector()
   {

   }
}

#endif