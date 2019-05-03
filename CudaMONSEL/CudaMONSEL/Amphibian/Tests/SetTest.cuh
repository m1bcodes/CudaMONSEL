#ifndef _SET_TEST_CUH_
#define _SET_TEST_CUH_

#include <cuda_runtime.h>

namespace SetTest
{
   class SetTest
   {
   public:
      __device__ SetTest();
      __device__ void TestIntBasic();
      __device__ void TestInt();
      __device__ void TestInt2();
      __device__ void TestString();
      __device__ void TestSetOfSetOfString();
   };

}

#endif
