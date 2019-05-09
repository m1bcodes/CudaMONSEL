#include "gov\nist\microanalysis\EPQTests\SphereTest.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace SphereTest
{
   static double SCALE = 1.0e-5;
   static int ITERATIONS = 1000;

   static PositionVecT makeCenter()
   {
      PositionVecT res({ (Math2::random() - 0.5) * SCALE, (Math2::random() - 0.5) * SCALE, (Math2::random() - 0.5) * SCALE });
      return res;
   }

   static PositionVecT makeNormal()
   {
      PositionVecT res({ Math2::random(), Math2::random(), Math2::random() });
      return Math2::normalize(res);
   }

   static void assertTrue(bool expr) {
      if (!expr) printf("false!");
   }

   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);
   }

   void testContains()
   {
      for (int i = 0; i < ITERATIONS; ++i) {
         auto center = makeCenter();
         double r = Math2::random() * SCALE + SCALE / 100.0;
         SphereT sphere(center.data(), r);
         {
            auto testPt = Math2::plus(center, Math2::multiply(0.99 * r * Math2::random(), makeNormal()));
            assertTrue(sphere.contains(testPt.data()));
         }
         {
            auto testPt = Math2::plus(center, Math2::multiply(r * (1.01 + Math2::random()), makeNormal()));
            assertTrue(!sphere.contains(testPt.data()));
         }
      }

      printf("SphereTest::testContains() completed.\n");
   }

   void testGetFirstIntersection() {
      for (int i = 0; i < ITERATIONS; ++i) {
         auto center = makeCenter();
         double r = Math2::random() * SCALE + SCALE / 100.0;
         SphereT sphere(center.data(), r);
         auto inside = Math2::plus(center, Math2::multiply(0.99 * r * Math2::random(), makeNormal()));
         assertTrue(sphere.contains(inside.data()));
         auto outside = Math2::plus(center, Math2::multiply(r * (1.01 + Math2::random()), makeNormal()));
         assertTrue(!sphere.contains(outside.data()));
         double t = sphere.getFirstIntersection(inside.data(), outside.data());
         double tp = sphere.getFirstIntersection(outside.data(), inside.data());
         assertTrue(t < 1.0);
         assertTrue(tp < 1.0);
         assertEquals(1.0, tp + t, 1.0e-6);
      }

      printf("SphereTest::testGetFirstIntersection() completed.\n");
   }
};