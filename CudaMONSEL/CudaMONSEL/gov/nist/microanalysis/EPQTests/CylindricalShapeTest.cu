#include "gov\nist\microanalysis\EPQTests\CylindricalShapeTest.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace CylindricalShapeTest
{
   VectorXd transform(const double pts[], double phi, double theta, double psi, const double offset[])
   {
      return Transform3D::translate(Transform3D::rotate(pts, phi, theta, psi).data(), offset, false);
   }

   static const double scale = (Math2::random() + 1.0e-4) * 10.0e-6;
   static const double radius = (Math2::random() + 1.0e-4) * 10.0e-6;
   static const double phi = Math2::random() * Math2::PI;
   static const double theta = Math2::random() * Math2::PI;
   static const double psi = Math2::random() * Math2::PI;
   static const double offset[] = { scale * Math2::random(), scale * Math2::random(), scale * Math2::random() };
   static const double end0[] = { -scale, 0.0, 0.0 };
   static const double end1[] = { scale, 0.0, 0.0 };

   static const VectorXd end0Vec = transform(end0, phi, theta, psi, offset);
   static const VectorXd end1Vec = transform(end1, phi, theta, psi, offset);
   static CylindricalShapeT shape(end0Vec.data(), end1Vec.data(), radius);

   double closestPtOnAxis(double pt[])
   {
      auto b = shape.getEnd0();
      auto ab = Math2::minus(shape.getEnd1(), b);
      VectorXd ptVec(pt, pt + 3);
      double t = Math2::dot(Math2::minus(ptVec, b), ab) / Math2::dot(ab, ab);
      return t;
   }

   bool isOnCylinder(double pt[])
   {
      double t = closestPtOnAxis(pt);
      if ((t >= 0) && (t <= 1)) {
         auto axisPt = Math2::plus(shape.getEnd0(), Math2::multiply(t, Math2::minus(shape.getEnd1(), shape.getEnd0())));
         VectorXd ptVec(pt, pt + 3);
         return ::abs(Math2::distance(ptVec, axisPt) - radius) < radius * 1.0e-6;
      }
      else
         return false;
   }

   bool isOnEndCap(double pt[])
   {
      double t = closestPtOnAxis(pt);
      VectorXd axisPt;
      if (::abs(t) < 1.0e-6)
         axisPt = shape.getEnd0();
      else if (::abs(t - 1.0) < 1.0e-6)
         axisPt = shape.getEnd1();
      else
         return false;
      VectorXd ptVec(pt, pt + 3);
      return axisPt.empty() ? false : Math2::distance(ptVec, axisPt) < radius;
   }

   static void assertTrue(bool expr)
   {
      if (!expr) printf("false!\n");
   }

   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);
   }

   void testZero()
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
   void testOne()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { -scale / 2.0, -radius / 2.0, 0.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnCylinder(Math2::pointBetween(parms0, parms1, t).data()));

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertEquals(1.0, tp + t, 1.0e-6);

      printf("CylindricalShapeTest::testOne() completed.\n");
   }

   /**
   * Test going through from one side to the other
   */
   void testTwo()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { -scale / 2.0, -radius / 2.0, -2.0 * radius };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnCylinder(Math2::pointBetween(parms0, parms1, t).data()));
      VectorXd pt = Math2::pointBetween(parms0, parms1, t);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(isOnCylinder(Math2::pointBetween(parms1, parms0, tp).data()));
      VectorXd ptp = Math2::pointBetween(parms1, parms0, tp);

      assertEquals(t, tp, 1.0e-6);

      assertEquals(Math2::distance(parms0, pt), Math2::distance(parms1, ptp), 1.0e-12);

      printf("CylindricalShapeTest::testTwo() completed.\n");
   }

   /**
   * Test going through the end caps
   */
   void testThree()
   {
      double parms0Vec[] = { -2.0 * scale, -radius, radius };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { 0.0, 0.0, 0.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms0, parms1, t).data()));
      auto pt = Math2::pointBetween(parms0, parms1, t);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms1, parms0, tp).data()));
      auto ptp = Math2::pointBetween(parms1, parms0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parms0, pt) + Math2::distance(parms1, ptp), Math2::distance(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testThree() completed.\n");
   }

   /**
   * Test going through the end caps
   */
   void testFour()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, 1.5 * radius };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { scale / 2.0, 0.0, 0.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms0, parms1, t).data()));
      VectorXd pt = Math2::pointBetween(parms0, parms1, t);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms1, parms0, tp).data()));
      VectorXd ptp = Math2::pointBetween(parms1, parms0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parms0, pt) + Math2::distance(parms1, ptp), Math2::distance(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testFour() completed.\n");
   }

   /**
   * Test parallel to axes
   */
   void testFive()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, radius / 2.0 };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { scale / 2.0, radius / 2.0, radius / 2.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms0, parms1, t).data()));
      auto pt = Math2::pointBetween(parms0, parms1, t);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms1, parms0, tp).data()));
      auto ptp = Math2::pointBetween(parms1, parms0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parms0, pt) + Math2::distance(parms1, ptp), Math2::distance(parms0, parms1), 1.0e-12);

      printf("CylindricalShapeTest::testFive() completed.\n");
   }

   /**
   * Test misses (parallel to axis)
   */
   void testSix()
   {
      double parms0Vec[] = { 2.0 * scale, 2.0 * radius, radius / 2.0 };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { -2.0 * scale, 2.0 * radius, radius / 2.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertEquals(t, INFINITY, 1.0e-6);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertEquals(tp, INFINITY, 1.0e-6);

      printf("CylindricalShapeTest::testSix() completed.\n");
   }

   /**
   * Test misses (parallel to axis)
   */
   void testSeven()
   {
      double parms0Vec[] = { 2.0 * scale, radius / 2.0, radius / 2.0 };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { 1.1 * scale, radius / 2.0, radius / 2.0 };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(t > 1.0);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(tp == INFINITY);

      printf("CylindricalShapeTest::testSeven() completed.\n");
   }

   /**
   * Test misses (not parallel to axis)
   */
   void testEight()
   {
      double parms0Vec[] = { -scale / 2.0, -radius / 2.0, 2.0 * radius };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { scale / 2.0, radius / 2.0, 1.1 * radius };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);
      double parms2Vec[] = { scale / 2.0, radius / 2.0, radius / 2.0 };
      VectorXd parm2 = transform(parms2Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(t > 1.0);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(tp > 1.0);

      double t2 = shape.getFirstIntersection(parms0.data(), parm2.data());
      assertTrue(isOnCylinder(Math2::pointBetween(parms0, parm2, t2).data()));

      double tp2 = shape.getFirstIntersection(parm2.data(), parms0.data());
      assertEquals(1.0, tp2 + t2, 1.0e-6);

      printf("CylindricalShapeTest::testEight() completed.\n");
   }

   /**
   * Test through both end cap and side (end0)
   */
   void testNine()
   {
      double parms0Vec[] = { -1.1 * scale, -radius / 10.0, 0.0 };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = { 0.0, 0.0, 1.1 * radius };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(isOnCylinder(Math2::pointBetween(parms1, parms0, tp).data()));

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(isOnEndCap(Math2::pointBetween(parms0, parms1, t).data()));

      assertTrue(1.0 + (tp + t) > 1.0e-6);

      printf("CylindricalShapeTest::testNine() completed.\n");
   }

   /**
   * Test through both end cap and side (end1)
   */
   void testTen()
   {
      double parms0Vec[] = {
         1.1 * scale,
         -radius / 10.0,
         0.0
      };
      VectorXd parms0 = transform(parms0Vec, phi, theta, psi, offset);
      double parms1Vec[] = {
         0.0,
         0.0,
         1.1 * radius
      };
      VectorXd parms1 = transform(parms1Vec, phi, theta, psi, offset);

      double t = shape.getFirstIntersection(parms0.data(), parms1.data());
      assertTrue(!isOnCylinder(Math2::pointBetween(parms0, parms1, t).data()));
      assertTrue(isOnEndCap(Math2::pointBetween(parms0, parms1, t).data()));

      double tp = shape.getFirstIntersection(parms1.data(), parms0.data());
      assertTrue(!isOnEndCap(Math2::pointBetween(parms1, parms0, tp).data()));
      assertTrue(isOnCylinder(Math2::pointBetween(parms1, parms0, tp).data()));

      assertTrue(1.0 - (tp + t) > 1.0e-6);

      printf("CylindricalShapeTest::testTen() completed.\n");
   }

   void testEleven()
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
         double r = 0.49 * SCALE * Math2::random();
         double th = Math2::random() * Math2::PI * 2.0;
         double inside[] = {
            1.9 * SCALE * (Math2::random() - 0.5),
            SCALE + ::cos(th) * r,
            SCALE + ::sin(th) * r
         };
         assertTrue(shape.contains(inside));
         th = Math2::random() * Math2::PI * 2.0;
         r = SCALE * Math2::random();
         double outside[] = {
            3.0 * SCALE * (Math2::random() - 0.5),
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

   void testTwelve()
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
      VectorXd saVec(sa, sa + 3);
      VectorXd sbVec(sb, sb + 3);
      auto pt = Math2::pointBetween(saVec, sbVec, t);
      assertEquals(::sqrt(Math2::sqr(pt[1]) + Math2::sqr(pt[2] - 1.0e-6)), 0.5e-6, 1.0e-12);

      printf("CylindricalShapeTest::testTwelve() completed.\n");
   }
}
