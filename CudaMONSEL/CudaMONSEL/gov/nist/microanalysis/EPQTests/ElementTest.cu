#include "ElementTest.cuh"
#include "..\EPQLibrary\Element.cuh"
#include "..\EPQLibrary\ToSI.cuh"

#include <stdio.h>

extern __device__ double fabs(double);

namespace ElementTest
{
   __device__ void assertEquals(double src, double target, double delta)
   {
      bool b = fabs(src - target) < delta;
      if (!b) {
         printf("values are different: %lf, %lf\n", src, target);
      }
   }

   __device__ void assertEquals(int src, int target)
   {
      if (src != target) {
         printf("values are different: %d, %d\n", src, target);
      }
   }

   __device__ void assertTrue(bool expr)
   {
      if (!expr) {
         printf("expr is not true\n");
      }
   }

   //__device__ void ElementTest::testOne()
   //{
   //   const Element elm = Element.byAtomicNumber(Element.elmTi);
   //   assertEquals(FromSI.eV(elm.meanIonizationPotential()), 247.24, 1.0);
   //   assertEquals(elm.getAtomicWeight(), 47.9, 1.0e-1);
   //   assertEquals(elm.getAtomicNumber(), 22);
   //   assertTrue(elm.compareTo(Element.byAtomicNumber(Element.elmH)) > 0);
   //   assertTrue(elm.compareTo(Element.byAtomicNumber(Element.elmFe)) < 0);
   //   assertTrue(elm.compareTo(Element.byAtomicNumber(Element.elmTi)) == 0);
   //   assertTrue(elm.equals(Element.byName("Ti")));
   //}
}
