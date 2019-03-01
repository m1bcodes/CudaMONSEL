#ifndef _MATERIAL_TEST_CUH_
#define _MATERIAL_TEST_CUH_

#include <cuda_runtime.h>

namespace MaterialTest
{
   class MaterialTest
   {
   public:
      __device__ MaterialTest();
      __device__ void testOne();
   };
}

#endif