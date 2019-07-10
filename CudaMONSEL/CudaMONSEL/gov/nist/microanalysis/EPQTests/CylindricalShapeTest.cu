#include "gov\nist\microanalysis\EPQTests\CylindricalShapeTest.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

#include "Amphibian\random.cuh"

namespace CylindricalShapeTest
{
   void transform3d(const double pts[], double phi, double theta, double psi, const double offset[], double res[])
   {
      double rotated[3];
      Transform3D::rotate3d(pts, phi, theta, psi, rotated);
      Transform3D::translate3d(rotated, offset, false, res);
   }

   static const double scale = (Random::random() + 1.0e-4) * 10.0e-6;
   static const double radius = (Random::random() + 1.0e-4) * 10.0e-6;
   static const double phi = Random::random() * Math2::PI;
   static const double theta = Random::random() * Math2::PI;
   static const double psi = Random::random() * Math2::PI;
   static const double offset[] = { scale * Random::random(), scale * Random::random(), scale * Random::random() };
   static const double end0[] = { -scale, 0.0, 0.0 };
   static const double end1[] = { scale, 0.0, 0.0 };

   double CylindricalShapeTest::closestPtOnAxis(const double pt[])
   {
      const double* b = shape.getEnd0();
      double ab[3];
      Math2::minus3d(shape.getEnd1(), b, ab);
      double m0[3];
      Math2::minus3d(pt, b, m0);
      return Math2::dot3d(m0, ab) / Math2::dot3d(ab, ab);
   }

   bool CylindricalShapeTest::isOnCylinder(const double pt[])
   {
      double t = closestPtOnAxis(pt);
      if ((t >= 0) && (t <= 1)) {
         double m0[3];
         Math2::minus3d(shape.getEnd1(), shape.getEnd0(), m0);
         double prd0[3];
         Math2::multiply3d(t, m0, prd0);
         double axisPt[3];
         Math2::plus3d(shape.getEnd0(), prd0, axisPt);
         return ::abs(Math2::distance3d(pt, axisPt) - radius) < radius * 1.0e-6;
      }
      else
         return false;
   }

   bool CylindricalShapeTest::isOnEndCap(const double pt[])
   {
      double t = closestPtOnAxis(pt);
      const double* axisPt;
      if (::abs(t) < 1.0e-6)
         axisPt = shape.getEnd0();
      else if (::abs(t - 1.0) < 1.0e-6)
         axisPt = shape.getEnd1();
      else
         return false;
      //return axisPt.empty() ? false : Math2::distance(ptVec, axisPt) < radius;
      return Math2::distance3d(pt, axisPt) < radius;
   }

   static void assertTrue(bool expr)
   {
      if (!expr) printf("false!\n");
   }

   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);
   }

   static double end0Transformed[3];
   static double end1Transformed[3];

   static double* transform3d2(const double pts[], double phi, double theta, double psi, const double offset[], int endidx)
   {
      if (endidx == 0) {
         transform3d(end0, phi, theta, psi, offset, end0Transformed);
         return end0Transformed;
      }
      else if (endidx == 1) {
         transform3d(end1, phi, theta, psi, offset, end1Transformed);
         return end1Transformed;
      }
      return nullptr;
   }

   CylindricalShapeTest::CylindricalShapeTest() :
      shape(transform3d2(end1, phi, theta, psi, offset, 1), transform3d2(end0, phi, theta, psi, offset, 0), radius)
   {
   }

   void CylindricalShapeTest::testZero()
   {
      double pts[] = { 0.5, 0.6, 0.7 };
      double pts2[] = { 1, 2, 3 };
      MatrixXd ret = Transform3D::rot(1, 2, 3);
      auto res0 = Transform3D::rotate(pts, 1, 2, 3);
      auto res = Transform3D::translate(pts, pts2, false);

      assertEquals(5.4030230587e-001, ::cos(1), 1e-6);
      assertEquals(8.4147098481e-001, ::sin(1), 1e-6);

      assertEquals(0.1038465651516682, ret[0][0], 1e-8);
      assertEquals(-0.42291857174254777, ret[0][1], 1e-8);
      assertEquals(-0.9001976297355174, ret[0][2], 1e-8);
      assertEquals(-0.864780102737098, ret[1][0], 1e-8);
      assertEquals(-0.4854784609636683, ret[1][1], 1e-8);
      assertEquals(0.12832006020245673, ret[1][2], 1e-8);
      assertEquals(-0.49129549643388193, ret[2][0], 1e-8);
      assertEquals(0.7651474012342926, ret[2][1], 1e-8);
      assertEquals(-0.4161468365471424, ret[2][2], 1e-8);

      assertEquals(-8.3196620128e-001, res0[0], 1e-8);
      assertEquals(-6.3385308581e-001, res0[1], 1e-8);
      assertEquals(-7.7862093059e-002, res0[2], 1e-8);

      assertEquals(1.5, res[0], 1e-8);
      assertEquals(2.6, res[1], 1e-8);
      assertEquals(3.7, res[2], 1e-8);

      printf("CylindricalShapeTest::testZero() completed.\n");
   }

   /**
   * Test going into and coming out of a side...
   */
   void CylindricalShapeTest::testOne()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { -scale / 2.0, -radius / 2.0, 0.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double ptbtw[3];
      Math2::pointBetween3d(parms0, parms1, t, ptbtw);
      assertTrue(isOnCylinder(ptbtw));

      double tp = shape.getFirstIntersection(parms1, parms0);
      assertEquals(1.0, tp + t, 1.0e-6);

      printf("CylindricalShapeTest::testOne() completed.\n");
   }

   /**
   * Test going through from one side to the other
   */
   void CylindricalShapeTest::testTwo()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { -scale / 2.0, -radius / 2.0, -2.0 * radius };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double pt[3];
      Math2::pointBetween3d(parms0, parms1, t, pt);
      assertTrue(isOnCylinder(pt));

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptp[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptp);
      assertTrue(isOnCylinder(ptp));

      assertEquals(t, tp, 1.0e-6);

      assertEquals(Math2::distance3d(parms0, pt), Math2::distance3d(parms1, ptp), 1.0e-12);

      printf("CylindricalShapeTest::testTwo() completed.\n");
   }

   /**
   * Test going through the end caps
   */
   void CylindricalShapeTest::testThree()
   {
      double parms0Vec[] = { -2.0 * scale, -radius, radius };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { 0.0, 0.0, 0.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double pt[3];
      Math2::pointBetween3d(parms0, parms1, t, pt);
      assertTrue(isOnEndCap(pt));

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptp[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptp);
      assertTrue(isOnEndCap(ptp));

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance3d(parms0, pt) + Math2::distance3d(parms1, ptp), Math2::distance3d(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testThree() completed.\n");
   }

   /**
   * Test going through the end caps
   */
   void CylindricalShapeTest::testFour()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, 1.5 * radius };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { scale / 2.0, 0.0, 0.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double pt[3];
      Math2::pointBetween3d(parms0, parms1, t, pt);
      assertTrue(isOnEndCap(pt));

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptp[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptp);
      assertTrue(isOnEndCap(ptp));

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance3d(parms0, pt) + Math2::distance3d(parms1, ptp), Math2::distance3d(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testFour() completed.\n");
   }

   /**
   * Test parallel to axes
   */
   void CylindricalShapeTest::testFive()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, radius / 2.0 };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { scale / 2.0, radius / 2.0, radius / 2.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double pt[3];
      Math2::pointBetween3d(parms0, parms1, t, pt);
      assertTrue(isOnEndCap(pt));

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptp[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptp);
      assertTrue(isOnEndCap(ptp));

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance3d(parms0, pt) + Math2::distance3d(parms1, ptp), Math2::distance3d(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testFive() completed.\n");
   }

   /**
   * Test misses (parallel to axis)
   */
   void CylindricalShapeTest::testSix()
   {
      double parms0Vec[] = { 2.0 * scale, 2.0 * radius, radius / 2.0 };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { -2.0 * scale, 2.0 * radius, radius / 2.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      assertEquals(t, INFINITY, 1.0e-6);

      double tp = shape.getFirstIntersection(parms1, parms0);
      assertEquals(tp, INFINITY, 1.0e-6);

      printf("CylindricalShapeTest::testSix() completed.\n");
   }

   /**
   * Test misses (parallel to axis)
   */
   void CylindricalShapeTest::testSeven()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, radius / 2.0 };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { 1.1 * scale, radius / 2.0, radius / 2.0 };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      assertTrue(t > 1.0);

      double tp = shape.getFirstIntersection(parms1, parms0);
      assertTrue(tp == INFINITY);

      printf("CylindricalShapeTest::testSeven() completed.\n");
   }

   /**
   * Test misses (not parallel to axis)
   */
   void CylindricalShapeTest::testEight()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { scale / 2.0, radius / 2.0, 1.1 * radius };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);
      double parms2Vec[] = { scale / 2.0, radius / 2.0, radius / 2.0 };
      double parms2[3];
      transform3d(parms2Vec, phi, theta, psi, offset, parms2);

      double t = shape.getFirstIntersection(parms0, parms1);
      assertTrue(t > 1.0);

      double tp = shape.getFirstIntersection(parms1, parms0);
      assertTrue(tp > 1.0);

      double t2 = shape.getFirstIntersection(parms0, parms2);
      double ptbtw[3];
      Math2::pointBetween3d(parms0, parms2, t2, ptbtw);
      assertTrue(isOnCylinder(ptbtw));

      double tp2 = shape.getFirstIntersection(parms2, parms0);
      assertEquals(1.0, tp2 + t2, 1.0e-6);

      printf("CylindricalShapeTest::testEight() completed.\n");
   }

   /**
   * Test through both end cap and side (end0)
   */
   void CylindricalShapeTest::testNine()
   {
      double parms0Vec[] = { -1.1 * scale, -radius / 10.0, 0.0 };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = { 0.0, 0.0, 1.1 * radius };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptbtw0[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptbtw0);
      assertTrue(isOnCylinder(ptbtw0));

      double t = shape.getFirstIntersection(parms0, parms1);
      double ptbtw1[3];
      Math2::pointBetween3d(parms0, parms1, t, ptbtw1);
      assertTrue(isOnEndCap(ptbtw1));

      assertTrue(1.0 + (tp + t) > 1.0e-6);

      printf("CylindricalShapeTest::testNine() completed.\n");
   }

   /**
   * Test through both end cap and side (end1)
   */
   void CylindricalShapeTest::testTen()
   {
      double parms0Vec[] = {
         1.1 * scale,
         -radius / 10.0,
         0.0
      };
      double parms0[3];
      transform3d(parms0Vec, phi, theta, psi, offset, parms0);
      double parms1Vec[] = {
         0.0,
         0.0,
         1.1 * radius
      };
      double parms1[3];
      transform3d(parms1Vec, phi, theta, psi, offset, parms1);

      double t = shape.getFirstIntersection(parms0, parms1);
      double ptbtw0[3];
      Math2::pointBetween3d(parms0, parms1, t, ptbtw0);
      assertTrue(!isOnCylinder(ptbtw0));
      double ptbtw1[3];
      Math2::pointBetween3d(parms0, parms1, t, ptbtw1);
      assertTrue(isOnEndCap(ptbtw1));

      double tp = shape.getFirstIntersection(parms1, parms0);
      double ptbtw2[3];
      Math2::pointBetween3d(parms1, parms0, tp, ptbtw2);
      assertTrue(!isOnEndCap(ptbtw2));
      assertTrue(isOnCylinder(ptbtw2));

      assertTrue(1.0 - (tp + t) > 1.0e-6);

      printf("CylindricalShapeTest::testTen() completed.\n");
   }

   void CylindricalShapeTest::testEleven()
   {
      double SCALE = 1.0e-5;
      int ITERATIONS = 1000;
      double p0[] = {
         -SCALE,
         SCALE,
         SCALE
      };
      double p1[] = {
         SCALE,
         SCALE,
         SCALE
      };

      CylindricalShapeT shape(p0, p1, 0.5 * SCALE);
      for (int i = 0; i < ITERATIONS; ++i) {
         double r = 0.49 * SCALE * Random::random();
         double th = Random::random() * Math2::PI * 2.0;
         double inside[] = {
            1.9 * SCALE * (Random::random() - 0.5),
            SCALE + ::cos(th) * r,
            SCALE + ::sin(th) * r
         };
         assertTrue(shape.contains(inside));
         th = Random::random() * Math2::PI * 2.0;
         r = SCALE * Random::random();
         double outside[] = {
            3.0 * SCALE * (Random::random() - 0.5),
            SCALE + ::cos(th) * (0.501 * SCALE + r),
            SCALE + ::sin(th) * (0.501 * SCALE + r)
         };
         assertTrue(!shape.contains(outside));
         double t = shape.getFirstIntersection(inside, outside);
         double tp = shape.getFirstIntersection(outside, inside);
         assertTrue(t < 1.0);
         assertTrue(tp < 1.0);
         assertEquals(1.0, t + tp, 1.0e-6);
      }

      printf("CylindricalShapeTest::testEleven() completed.\n");
   }

   void CylindricalShapeTest::testTwelve()
   {
      double pt0[] = {
         1.0e-6,
         0.0,
         1.0e-6
      };
      double pt1[] = {
         -1.0e-6,
         0.0,
         1.0e-6
      };
      CylindricalShapeT shape(pt0, pt1, 0.5e-6);
      double sa[] = {
         -1.413972850134937E-7,
         -1.5600411637508016E-7,
         1.4006819741698632E-6
      };
      assertTrue(shape.contains(sa));
      double sb[] = {
         -8.248126103570508E-9,
         -2.5333627600912425E-7,
         7.734838104262905E-7
      };
      assertTrue(shape.contains(sb));
      double t = shape.getFirstIntersection(sa, sb);
      assertTrue(t != INFINITY);
      assertTrue(t > 1.0);
      double pt[3];
      Math2::pointBetween3d(sa, sb, t, pt);
      assertEquals(::sqrt(Math2::sqr(pt[1]) + Math2::sqr(pt[2] - 1.0e-6)), 0.5e-6, 1.0e-12);

      printf("CylindricalShapeTest::testTwelve() completed.\n");
   }
}
