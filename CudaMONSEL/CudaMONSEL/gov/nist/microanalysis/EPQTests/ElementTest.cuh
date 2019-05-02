#ifndef _ELEMENT_TEST_CUH_
#define _ELEMENT_TEST_CUH_

namespace ElementTest
{
   class ElementTest
   {
   public:
      ElementTest();
      void testZero();
      void testOne();
   };

   void assertEquals(double, double, double);
   void assertEquals(int src, int target);
   void assertTrue(bool);
}

#endif
