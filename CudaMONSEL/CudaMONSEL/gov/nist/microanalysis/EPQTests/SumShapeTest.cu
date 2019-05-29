#include "gov\nist\microanalysis\EPQTests\SumShapeTest.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\BasicMaterialModel.cuh"

//#include "gov\nist\microanalysis\NISTMonte\TrajectoryVRML.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace SumShapeTest
{
   static const ElementT* mat1Elements[] = {
      &Element::C,
      &Element::O,
      &Element::Cl
   };
   static const double mat1MassFracs[] = {
      0.7,
      0.28,
      0.02
   };
   //static const CompositionT mat1Comp(mat1Elements, 3, mat1MassFracs, 3, "Resin");
   static const MaterialT mat1(mat1Elements, 3, mat1MassFracs, 3, ToSI::gPerCC(1.14), "Resin");

   static const ElementT* mat2Element[] = {
      &Element::C,
      &Element::O,
      &Element::N,
      &Element::H,
      &Element::S,
      &Element::P,
      &Element::Os
   };
   static const double mat2MassFracs[] = {
      0.4962,
      0.2988,
      0.0784,
      0.0773,
      0.006,
      0.0232,
      0.02
   };
   //static const CompositionT mat2Comp(mat2Element, 7, mat2MassFracs, 7, "Inner");
   static const MaterialT mat2(mat2Element, 7, mat2MassFracs, 7, ToSI::gPerCC(1.11), "Inner");

   static const double ChamberRadius = 0.1;

   // The chamber is a Sphere centered at the origin with 0.10 m radius
   static const double center[] = {
      0.0,
      0.0,
      0.0
   };
   static SphereT sphere(center, ChamberRadius);
   static NullMaterialScatterModelT NULL_MSM;
   static RegionT mChamber(nullptr, &NULL_MSM, &sphere);

   static const double beamCenter[] = {
      0,
      0,
      -0.05
   };
   static const double beamEnergy = ToSI::keV(20.0);
   static const double beamDia = 1.0e-8;
   static const GaussianBeamT beam(beamDia, beamEnergy, beamCenter);

   static MonteCarloSST mMonte(&beam, &mChamber, beam.createElectron());

   static PlaneT pl(Math2::MINUS_Z_AXIS, 3, Math2::ORIGIN_3D, 3);
   static PlaneT* pls[] = { &pl };
   static MultiPlaneShape::MultiPlaneShape blk(pls, 1);

   static const double radius = 0.5e-6;
   static const double length = 1.0e-6;
   static const double end0[] = {
      -length,
      0.0,
      2.0 * radius
   };
   static const double end1[] = {
      length,
      0.0,
      2.0 * radius
   };

   static SphereT mCap0Outer(end0, radius);
   static SphereT mCap1Outer(end1, radius);
   static CylindricalShapeT mCylOuter(end0, end1, radius);

   static ShapeT* const shapes[] = {
      &mCylOuter,
      &mCap0Outer,
      &mCap1Outer
   };

   static SumShapeT mPillOuter(shapes, 3);

   static BasicMaterialModelT bmm1(mat1);
   static RegionT r1(&mChamber, &bmm1, &blk);

   static BasicMaterialModelT bmm2(mat2);
   static RegionT r2(&r1, &bmm2, &mPillOuter);

   static void assertTrue(bool expr)
   {
      if (!expr) printf("false!\n");
   }

   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);;
   }

   static VectorXd pointInside()
   {
      VectorXd res(3);
      switch (Math2::randomInt(3)) {
      case 0: {
         const double r = radius * Math2::random();
         const double th = Math2::random() * Math2::PI;
         const double phi = 2.0 * Math2::random() * Math2::PI;
         res[0] = -length + r * ::cos(phi) * ::sin(th);
         res[1] = r * ::sin(phi) * ::sin(th);
         res[2] = 2.0 * radius + r * ::cos(th);
         assertTrue(mCap0Outer.contains(res.data()));
      }
              break;
      case 1: {
         const double phi = 2.0 * Math2::random() * Math2::PI;
         const double r = radius * Math2::random();
         res[0] = 2.0 * length * (Math2::random() - 0.5);
         res[1] = r * ::cos(phi);
         res[2] = 2.0 * radius + ::sin(phi) * r;
         assertTrue(mCylOuter.contains(res.data()));
      }
              break;
      case 2: {
         const double r = radius * Math2::random();
         const double th = Math2::random() * Math2::PI;
         const double phi = 2.0 * Math2::random() * Math2::PI;
         res[0] = length + r * ::cos(phi) * ::sin(th);
         res[1] = r * ::sin(phi) * ::sin(th);
         res[2] = 2.0 * radius + r * ::cos(th);
         assertTrue(mCap1Outer.contains(res.data()));
      }
              break;
      }
      return res;
   }

   static VectorXd pointOutside()
   {
      VectorXd res(3);
      const double phi = 2.0 * Math2::random() * Math2::PI;
      const double r = radius * (1.00001 + 0.9 * Math2::random());
      res[0] = 3.0 * length * (Math2::random() - 0.5);
      res[1] = r * ::cos(phi);
      res[2] = 2.0 * radius + ::sin(phi) * r;
      return res;
   }

   void testGetFirstIntersection()
   {
      for (int i = 0; i < 1000; ++i) {
         auto inside1 = pointInside();
         assertTrue(mPillOuter.contains(inside1.data()));
         auto inside2 = pointInside();
         assertTrue(mPillOuter.contains(inside2.data()));
         auto outside = pointOutside();
         assertTrue(!mPillOuter.contains(outside.data()));
         {
            const double t = mPillOuter.getFirstIntersection(inside1.data(), inside2.data());
            assertTrue(t != INFINITY);
            assertTrue(t > 1.0);
         }
         {
            //TrajectoryVRML vrml = null;
            //PrintWriter wr = null;
            //if (false)
            //   try {
            //   wr = new PrintWriter(new File("C:\\Documents and Settings\\Nicholas\\Desktop\\dump.wrl"));
            //   vrml = new TrajectoryVRML(mMonte, wr);
            //   vrml.renderSample();
            //   wr.println("# full line");
            //   vrml.drawLine(inside1, outside);
            //}
            //catch (const Exception e) {
            //   e.printStackTrace();
            //}
            assertTrue(mPillOuter.contains(inside1.data()));
            assertTrue(!mPillOuter.contains(outside.data()));
            const double t = mPillOuter.getFirstIntersection(inside1.data(), outside.data());
            //if (vrml != null) {
            //   if (wr != null)
            //      wr.println("# inside to outside");
            //   vrml.drawLine(inside1, Math2::pointBetween(inside1, outside, t));
            //}
            assertTrue(mPillOuter.contains(Math2::pointBetween(inside1, outside, 0.99 * t).data()));
            assertTrue(!mPillOuter.contains(Math2::pointBetween(inside1, outside, 1.01 * t).data()));
            const double tp = mPillOuter.getFirstIntersection(outside.data(), inside1.data());
            //if (vrml != null) {
            //   wr.println("# outside to inside");
            //   vrml.drawLine(outside, Math2::pointBetween(outside, inside1, tp));
            //   wr.flush();
            //}
            assertTrue(!mPillOuter.contains(Math2::pointBetween(outside, inside1, 0.99 * tp).data()));
            assertTrue(mPillOuter.contains(Math2::pointBetween(outside, inside1, 1.01 * tp).data()));
            assertTrue(t < 1.0);
            assertTrue(tp < 1.0);
            assertEquals(1.0, t + tp, 1.0e-6);
         }
         {
            const double t = mPillOuter.getFirstIntersection(inside2.data(), outside.data());
            const double tp = mPillOuter.getFirstIntersection(outside.data(), inside2.data());
            assertTrue(t < 1.0);
            assertTrue(tp < 1.0);
            assertEquals(1.0, t + tp, 1.0e-6);
         }
      }

      printf("SumShapeTest::testGetFirstIntersection() completed.\n");
   }

   void testAll()
   {
      mMonte.runMultipleTrajectories(1000);

      printf("SumShapeTest::testAll() completed.\n");
   }
}
