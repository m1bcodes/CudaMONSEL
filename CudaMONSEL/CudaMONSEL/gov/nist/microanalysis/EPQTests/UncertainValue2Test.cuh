#ifndef _UNCERTAIN_VALUES_2_TEST_CUH_
#define _UNCERTAIN_VALUES_2_TEST_CUH_

#include "..\Utility\UncertainValue2.cuh"

namespace UncertainValue2Test
{
   class UncertainValue2Test
   {
   public:
      __device__ UncertainValue2Test();
      __device__ void testSpecialValues();
      __device__ void testA();
      __device__ void testB();
      __device__ void testC();
      __device__ void testAB();
      __device__ void testAdd1();
      __device__ void testAdd2();
      __device__ void testAdd3();
      __device__ void testMultiply();
      __device__ void testDivide();
      __device__ void testFunctions();
   };

   __device__ UncertainValue2::UncertainValue2 makeA2a();
   __device__ UncertainValue2::UncertainValue2 makeB2a();
   __device__ UncertainValue2::UncertainValue2 makeC2a();
   __device__ void assertEquals(UncertainValue2::UncertainValue2 uv, UncertainValue2::UncertainValue2 uv2, double delta);
}
#endif