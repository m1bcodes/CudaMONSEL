#include "gov\nist\microanalysis\EPQTests\SphereTest.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace SphereTest
{
   static double SCALE = 1.0e-5;
   static int ITERATIONS = 1000;

   static void makeCenter(double res[])
   {
      res[0] = (Random::random() - 0.5) * SCALE;
      res[1] = (Random::random() - 0.5) * SCALE;
      res[2] = (Random::random() - 0.5) * SCALE;
   }

   static void makeNormal(double res[])
   {
      double norm[] = { Random::random(), Random::random(), Random::random() };
      Math2::normalize3d(norm, res);
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
         double center[3];
         makeCenter(center);
         double r = Random::random() * SCALE + SCALE / 100.0;
         SphereT sphere(center, r);
         {
            double norm[3];
            makeNormal(norm);
            double mult[3];
            Math2::multiply3d(0.99 * r * Random::random(), norm, mult);
            double testPt[3];
            Math2::plus3d(center, mult, testPt);
            assertTrue(sphere.contains(testPt));
         }
         {
            double norm[3];
            makeNormal(norm);
            double mult[3];
            Math2::multiply3d(r * (1.01 + Random::random()), norm, mult);
            double testPt[3];
            Math2::plus3d(center, mult, testPt);
            assertTrue(!sphere.contains(testPt));
         }
      }

      printf("SphereTest::testContains() completed.\n");
   }

   void testGetFirstIntersection() {
      double center[3];
      double norm[3];
      double mult[3];
      for (int i = 0; i < ITERATIONS; ++i) {
         makeCenter(center);
         double r = Random::random() * SCALE + SCALE / 100.0;
         SphereT sphere(center, r);
         makeNormal(norm);
         Math2::multiply3d(0.99 * r * Random::random(), norm, mult);
         double inside[3];
         Math2::plus3d(center, mult, inside);
         assertTrue(sphere.contains(inside));
         makeNormal(norm);
         Math2::multiply3d(r * (1.01 + Random::random()), norm, mult);
         double outside[3];
         Math2::plus3d(center, mult, outside);
         assertTrue(!sphere.contains(outside));
         double t = sphere.getFirstIntersection(inside, outside);
         double tp = sphere.getFirstIntersection(outside, inside);
         assertTrue(t < 1.0);
         assertTrue(tp < 1.0);
         assertEquals(1.0, tp + t, 1.0e-6);
      }

      printf("SphereTest::testGetFirstIntersection() completed.\n");
   }
};