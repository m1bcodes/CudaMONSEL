#ifndef _UNCERTAIN_VALUES_2_TEST_CUH_
#define _UNCERTAIN_VALUES_2_TEST_CUH_

#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"

namespace UncertainValue2Test
{
   class UncertainValue2Test
   {
   public:
      UncertainValue2Test();

      void testSpecialValues();
      void testA();
      void testB();
      void testC();
      void testAB();
      void testAdd1();
      void testAdd2();
      void testAdd3();
      void testMultiply();
      void testDivide();
      void testFunctions();
   };

   UncertainValue2::UncertainValue2 makeA2a();
   UncertainValue2::UncertainValue2 makeB2a();
   UncertainValue2::UncertainValue2 makeC2a();

   void assertEquals(UncertainValue2::UncertainValue2 uv, UncertainValue2::UncertainValue2 uv2, double delta);
}
#endif