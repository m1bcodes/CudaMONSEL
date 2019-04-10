#include "ElementTest.cuh"
#include "..\EPQLibrary\Element.cuh"
#include "..\EPQLibrary\FromSI.cuh"

#include <stdio.h>
#include <math.h>

namespace ElementTest
{
   void assertEquals(double src, double target, double delta)
   {
      bool b = fabs(src - target) < delta;
      if (!b) {
         printf("ElementTest::assertEquals: values are different: %lf, %lf\n", src, target);
      }
   }

   void assertEquals(int src, int target)
   {
      if (src != target) {
         printf("ElementTest::assertEquals: values are different: %d, %d\n", src, target);
      }
   }

   void assertTrue(bool expr)
   {
      if (!expr) {
         printf("ElementTest::assertTrue: expr is not true\n");
      }
   }

   ElementTest::ElementTest()
   {
   }

   void ElementTest::testOne()
   {
      Element::Element elm = Element::byAtomicNumber(Element::elmTi);
      //assertEquals(FromSI::eV(elm.meanIonizationPotential()), 247.24, 1.0);
      assertEquals(elm.getAtomicWeight(), 47.9, 1.0e-1);
      //assertEquals(elm.getIonizationEnergy(), 6.8281, 1.0e-1);
      //assertEquals(elm.getAtomicNumber(), 22);
      //assertTrue(elm.compareTo(Element::byAtomicNumber(Element::elmH)) > 0);
      //assertTrue(elm.compareTo(Element::byAtomicNumber(Element::elmFe)) < 0);
      //assertTrue(elm.compareTo(Element::byAtomicNumber(Element::elmTi)) == 0);
      //assertTrue(elm.equals(Element::byName("Ti")));

      printf("%s completed.\n", "ElementTest::testOne()");
   }
}
