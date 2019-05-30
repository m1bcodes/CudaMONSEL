// file: gov\nist\microanalysis\EPQTests\CylindricalShapeTest.cuh

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"

#ifndef _CYLINDRICAL_SHAPE_TEST_CUH_
#define _CYLINDRICAL_SHAPE_TEST_CUH_

namespace CylindricalShapeTest
{
   class CylindricalShapeTest
   {
   public:
      CylindricalShapeTest();

      void testZero();
      void testOne();
      void testTwo();
      void testThree();
      void testFour();
      void testFive();
      void testSix();
      void testSeven();
      void testEight();
      void testNine();
      void testTen();
      void testEleven();
      void testTwelve();

   private:
      double closestPtOnAxis(const double pt[]);
      bool isOnCylinder(const double pt[]);
      bool isOnEndCap(const double pt[]);

      CylindricalShapeT shape;
   };
}

#endif