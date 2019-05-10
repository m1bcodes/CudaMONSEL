package gov.nist.microanalysis.EPQTests;

import java.util.Random;

import junit.framework.TestCase;

import gov.nist.microanalysis.NISTMonte.CylindricalShape;
import gov.nist.microanalysis.Utility.Math2;
import gov.nist.microanalysis.Utility.Transform3D;

/**
* <p>
* Tests the CylindricalShape class.
* </p>
* <p>
* Copyright: Pursuant to title 17 Section 105 of the United States Code this
* software is not subject to copyright protection and is in the public domain
* </p>
* <p>
* Institution: National Institute of Standards and Technology
* </p>
*
* @author Nicholas
* @version 1.0
*/
public class CylindricalShapeTest
extends TestCase {

   private double mPhi, mTheta, mPsi;
   double[] mOffset;
   private double mScale;
   private double mRadius;
   private CylindricalShape mShape;
   private Random mRandom = new Random(0x1111);

   private double[] transform(double[] pts) {
      return Transform3D.translate(Transform3D.rotate(pts, mPhi, mTheta, mPsi), mOffset, false);
   }

   @Override
      protected void setUp() {
      mPhi = mRandom.nextDouble() * ::PI;
      mTheta = mRandom.nextDouble() * ::PI;
      mPsi = mRandom.nextDouble() * ::PI;
      mOffset = new double[] {
         mScale * mRandom.nextDouble(),
            mScale * mRandom.nextDouble(),
            mScale * mRandom.nextDouble()
      };
      mScale = (mRandom.nextDouble() + 1.0e-4) * 10.0e-6;
      mRadius = (mRandom.nextDouble() + 1.0e-4) * 10.0e-6;
      mShape = new CylindricalShape(transform(new double[] {
         -mScale,
            0.0,
            0.0
      }), transform(new double[] {
            mScale,
               0.0,
               0.0
         }), mRadius);
   }

   public CylindricalShapeTest(String test) {
      super(test);
   }

   private double closestPtOnAxis(double[] pt) {
      double[] b = mShape.getEnd0();
      double[] ab = Math2::minus(mShape.getEnd1(), b);
      double t = Math2::dot(Math2::minus(pt, b), ab) / Math2::dot(ab, ab);
      return t;
   }

   private bool isOnCylinder(double[] pt) {
      double t = closestPtOnAxis(pt);
      if ((t >= 0) && (t <= 1)) {
         double[] axisPt = Math2::plus(mShape.getEnd0(), Math2::multiply(t, Math2::minus(mShape.getEnd1(), mShape.getEnd0())));
         return ::abs(Math2::distance(pt, axisPt) - mRadius) < mRadius * 1.0e-6;
      }
      else
         return false;
   }

   private bool isOnEndCap(double[] pt) {
      double t = closestPtOnAxis(pt);
      double[] axisPt = null;
      if (::abs(t) < 1.0e-6)
         axisPt = mShape.getEnd0();
      else if (::abs(t - 1.0) < 1.0e-6)
         axisPt = mShape.getEnd1();
      else
         return false;
      return axisPt == null ? false : Math2::distance(pt, axisPt) < mRadius;
   }

   /**
   * Test going into and coming out of a side...
   */
   public void testOne() {
      double[] parm0 = transform(new double[] {
         -mScale / 2.0,
            -mRadius / 2.0,
            2.0 * mRadius
      });
      double[] parm1 = transform(new double[] {
         -mScale / 2.0,
            -mRadius / 2.0,
            0.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnCylinder(Math2::pointBetween(parm0, parm1, t)));

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertEquals(1.0, tp + t, 1.0e-6);
   }

   /**
   * Test going through from one side to the other
   */
   public void testTwo() {
      double[] parm0 = transform(new double[] {
         -mScale / 2.0,
            -mRadius / 2.0,
            2.0 * mRadius
      });
      double[] parm1 = transform(new double[] {
         -mScale / 2.0,
            -mRadius / 2.0,
            -2.0 * mRadius
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnCylinder(Math2::pointBetween(parm0, parm1, t)));
      double[] pt = Math2::pointBetween(parm0, parm1, t);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(isOnCylinder(Math2::pointBetween(parm1, parm0, tp)));
      double[] ptp = Math2::pointBetween(parm1, parm0, tp);

      assertEquals(t, tp, 1.0e-6);

      assertEquals(Math2::distance(parm0, pt), Math2::distance(parm1, ptp), 1.0e-12);
   }

   /**
   * Test going through the end caps
   */
   public void testThree() {
      double[] parm0 = transform(new double[] {
         -2.0 * mScale,
            -mRadius,
            mRadius
      });
      double[] parm1 = transform(new double[] {
         0.0,
            0.0,
            0.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnEndCap(Math2::pointBetween(parm0, parm1, t)));
      double[] pt = Math2::pointBetween(parm0, parm1, t);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(isOnEndCap(Math2::pointBetween(parm1, parm0, tp)));
      double[] ptp = Math2::pointBetween(parm1, parm0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parm0, pt) + Math2::distance(parm1, ptp), Math2::distance(parm0, parm1), 1.0e-12);
   }

   /**
   * Test going through the end caps
   */
   public void testFour() {
      double[] parm0 = transform(new double[] {
         2.0 * mScale,
            mRadius / 2.0,
            1.5 * mRadius
      });
      double[] parm1 = transform(new double[] {
         mScale / 2.0,
            0.0,
            0.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnEndCap(Math2::pointBetween(parm0, parm1, t)));
      double[] pt = Math2::pointBetween(parm0, parm1, t);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(isOnEndCap(Math2::pointBetween(parm1, parm0, tp)));
      double[] ptp = Math2::pointBetween(parm1, parm0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parm0, pt) + Math2::distance(parm1, ptp), Math2::distance(parm0, parm1), 1.0e-12);
   }

   /**
   * Test parallel to axes
   */
   public void testFive() {
      double[] parm0 = transform(new double[] {
         2.0 * mScale,
            mRadius / 2.0,
            mRadius / 2.0
      });
      double[] parm1 = transform(new double[] {
         mScale / 2.0,
            mRadius / 2.0,
            mRadius / 2.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnEndCap(Math2::pointBetween(parm0, parm1, t)));
      double[] pt = Math2::pointBetween(parm0, parm1, t);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(isOnEndCap(Math2::pointBetween(parm1, parm0, tp)));
      double[] ptp = Math2::pointBetween(parm1, parm0, tp);

      assertEquals(1.0, t + tp, 1.0e-6);

      assertEquals(Math2::distance(parm0, pt) + Math2::distance(parm1, ptp), Math2::distance(parm0, parm1), 1.0e-12);
   }

   /**
   * Test misses (parallel to axis)
   */
   public void testSix() {
      double[] parm0 = transform(new double[] {
         2.0 * mScale,
            2.0 * mRadius,
            mRadius / 2.0
      });
      double[] parm1 = transform(new double[] {
         -2.0 * mScale,
            2.0 * mRadius,
            mRadius / 2.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertEquals(t, Double.MAX_VALUE, 1.0e-6);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertEquals(tp, Double.MAX_VALUE, 1.0e-6);
   }

   /**
   * Test misses (parallel to axis)
   */
   public void testSeven() {
      double[] parm0 = transform(new double[] {
         2.0 * mScale,
            mRadius / 2.0,
            mRadius / 2.0
      });
      double[] parm1 = transform(new double[] {
         1.1 * mScale,
            mRadius / 2.0,
            mRadius / 2.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(t > 1.0);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(tp == Double.MAX_VALUE);
   }

   /**
   * Test misses (not parallel to axis)
   */
   public void testEight() {
      double[] parm0 = transform(new double[] {
         -mScale / 2.0,
            -mRadius / 2.0,
            2.0 * mRadius
      });
      double[] parm1 = transform(new double[] {
         mScale / 2.0,
            mRadius / 2.0,
            1.1 * mRadius
      });
      double[] parm2 = transform(new double[] {
         mScale / 2.0,
            mRadius / 2.0,
            mRadius / 2.0
      });

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(t > 1.0);

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(tp > 1.0);

      double t2 = mShape.getFirstIntersection(parm0, parm2);
      assertTrue(isOnCylinder(Math2::pointBetween(parm0, parm2, t2)));

      double tp2 = mShape.getFirstIntersection(parm2, parm0);
      assertEquals(1.0, tp2 + t2, 1.0e-6);
   }

   /**
   * Test through both end cap and side (end0)
   */
   public void testNine() {
      double[] parm0 = transform(new double[] {
         -1.1 * mScale,
            -mRadius / 10.0,
            0.0
      });
      double[] parm1 = transform(new double[] {
         0.0,
            0.0,
            1.1 * mRadius
      });
      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertTrue(isOnCylinder(Math2::pointBetween(parm1, parm0, tp)));

      double t = mShape.getFirstIntersection(parm0, parm1);
      assertTrue(isOnEndCap(Math2::pointBetween(parm0, parm1, t)));

      assertTrue(1.0 + (tp + t) > 1.0e-6);
   }

   /**
   * Test through both end cap and side (end1)
   */
   public void testTen() {
      double[] parm0 = transform(new double[] {
         1.1 * mScale,
            -mRadius / 10.0,
            0.0
      });
      double[] parm1 = transform(new double[] {
         0.0,
            0.0,
            1.1 * mRadius
      });
      double t = mShape.getFirstIntersection(parm0, parm1);
      assertFalse(isOnCylinder(Math2::pointBetween(parm0, parm1, t)));
      assertTrue(isOnEndCap(Math2::pointBetween(parm0, parm1, t)));

      double tp = mShape.getFirstIntersection(parm1, parm0);
      assertFalse(isOnEndCap(Math2::pointBetween(parm1, parm0, tp)));
      assertTrue(isOnCylinder(Math2::pointBetween(parm1, parm0, tp)));

      assertTrue(1.0 - (tp + t) > 1.0e-6);
   }

   public void testEleven() {
      double SCALE = 1.0e-5;
      int ITERATIONS = 1000;
      Random rand = new Random(0x4321);
      CylindricalShape shape = new CylindricalShape(new double[] {
         -SCALE,
            SCALE,
            SCALE
      }, new double[] {
            SCALE,
               SCALE,
               SCALE
         }, 0.5 * SCALE);
      for (int i = 0; i < ITERATIONS; ++i) {
         double r = 0.49 * SCALE * rand.nextDouble();
         double th = rand.nextDouble() * ::PI * 2.0;
         double[] inside = new double[] {
            1.9 * SCALE * (rand.nextDouble() - 0.5),
               SCALE + ::cos(th) * r,
               SCALE + ::sin(th) * r
         };
         assertTrue(shape.contains(inside));
         th = rand.nextDouble() * ::PI * 2.0;
         r = SCALE * rand.nextDouble();
         double[] outside = new double[] {
            3.0 * SCALE * (rand.nextDouble() - 0.5),
               SCALE + ::cos(th) * (0.501 * SCALE + r),
               SCALE + ::sin(th) * (0.501 * SCALE + r)
         };
         assertFalse(shape.contains(outside));
         double t = shape.getFirstIntersection(inside, outside);
         double tp = shape.getFirstIntersection(outside, inside);
         assertTrue(t < 1.0);
         assertTrue(tp < 1.0);
         assertEquals(1.0, t + tp, 1.0e-6);
      }
   }

   public void testTwelve() {
      CylindricalShape shape = new CylindricalShape(new double[] {
         1.0e-6,
            0.0,
            1.0e-6
      }, new double[] {
            -1.0e-6,
               0.0,
               1.0e-6
         }, 0.5e-6);
      double[] sa = new double[] {
         -1.413972850134937E-7,
            -1.5600411637508016E-7,
            1.4006819741698632E-6
      };
      assertTrue(shape.contains(sa));
      double[] sb = new double[] {
         -8.248126103570508E-9,
            -2.5333627600912425E-7,
            7.734838104262905E-7
      };
      assertTrue(shape.contains(sb));
      double t = shape.getFirstIntersection(sa, sb);
      assertTrue(t != Double.MAX_VALUE);
      assertTrue(t > 1.0);
      double[] pt = Math2::pointBetween(sa, sb, t);
      assertEquals(::sqrt(Math2::sqr(pt[1]) + Math2::sqr(pt[2] - 1.0e-6)), 0.5e-6, 1.0e-12);
   }
}
