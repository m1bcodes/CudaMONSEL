#ifndef _SET_TEST_CUH_
#define _SET_TEST_CUH_

#include <cuda_runtime.h>

namespace SetTest
{
   class SetTest
   {
   public:
      __host__ __device__ SetTest();
      __host__ __device__ void testIntBasic();
      __host__ __device__ void testInt();
      __host__ __device__ void testInt2();
      __host__ __device__ void testString();
      //__device__ void TestSetOfSetOfString();
   };
}

#endif