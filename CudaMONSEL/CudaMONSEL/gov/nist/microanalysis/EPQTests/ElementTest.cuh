#ifndef _ELEMENT_TEST_CUH_
#define _ELEMENT_TEST_CUH_

#include <cuda_runtime.h>

namespace ElementTest
{
   //class ElementTest
   //{
   //public:
   //   __device__ void testOne();
   //};

   __device__ void assertEquals(double, double, double);
   __device__ void assertEquals(int src, int target);
   __device__ void assertTrue(bool);
}

#endif
