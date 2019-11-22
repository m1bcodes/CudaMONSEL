#include "gov\nist\microanalysis\EPQTests\SumShapeTest.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

//#include "gov\nist\microanalysis\NISTMonte\TrajectoryVRML.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\random.cuh"

namespace SumShapeTest
{
   static const ElementT* mat1Elements[] = {
      &Element::C,
      &Element::O,
      &Element::Cl
   };
   static const ::Composition::data_type mat1MassFracs[] = {
      0.7,
      0.28,
      0.02
   };
   //static const CompositionT mat1Comp(mat1Elements, 3, mat1MassFracs, 3, "Resin");
   //static const MaterialT mat1(mat1Elements, 3, mat1MassFracs, 3, ToSI::gPerCC(1.14), "Resin");

   static const ElementT* mat2Element[] = {
      &Element::C,
      &Element::O,
      &Element::N,
      &Element::H,
      &Element::S,
      &Element::P,
      &Element::Os
   };
   static const ::Composition::data_type mat2MassFracs[] = {
      0.4962,
      0.2988,
      0.0784,
      0.0773,
      0.006,
      0.0232,
      0.02
   };
   //static const CompositionT mat2Comp(mat2Element, 7, mat2MassFracs, 7, "Inner");
   //static const MaterialT mat2(mat2Element, 7, mat2MassFracs, 7, ToSI::gPerCC(1.11), "Inner");

   static const double ChamberRadius = 0.1;

   // The chamber is a Sphere centered at the origin with 0.10 m radius
   static const double center[] = {
      0.0,
      0.0,
      0.0
   };
   //static SphereT sphere(center, ChamberRadius);
   //static NullMaterialScatterModelT NULL_MSM;
   //static RegionT mChamber(nullptr, &NULL_MSM, &sphere);

   static const double beamCenter[] = {
      0,
      0,
      -0.05
   };
   static const double beamEnergy = ToSI::keV(20.0);
   static const double beamDia = 1.0e-8;
   //static const GaussianBeamT beam(beamDia, beamEnergy, beamCenter);

   //static MonteCarloSST mMonte(&beam, &mChamber, beam.createElectron());

   //static PlaneT pl(Math2::MINUS_Z_AXIS, Math2::ORIGIN_3D);
   //static PlaneT* pls[] = { &pl };
   //static MultiPlaneShapeT blk(pls, 1);

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

   //static SphereT mCap0Outer(end0, radius);
   //static SphereT mCap1Outer(end1, radius);
   //static CylindricalShapeT mCylOuter(end0, end1, radius);

   //static ShapeT* const shapes[] = {
   //   &mCylOuter,
   //   &mCap0Outer,
   //   &mCap1Outer
   //};

   //static SumShapeT mPillOuter(shapes, 3);

   //static BasicMaterialModelT bmm1(mat1);
   //static RegionT r1(&mChamber, &bmm1, &blk);

   //static BasicMaterialModelT bmm2(mat2);
   //static RegionT r2(&r1, &bmm2, &mPillOuter);

   SumShapeTest::SumShapeTest() :
      mat1(mat1Elements, 3, mat1MassFracs, 3, ToSI::gPerCC(1.14), "Resin"),
      mat2(mat2Element, 7, mat2MassFracs, 7, ToSI::gPerCC(1.11), "Inner"),
      sphere(center, ChamberRadius),
      mChamber(nullptr, &NULL_MSM, &sphere),
      beam(beamDia, beamEnergy, beamCenter),
      //mMonte(&beam, &mChamber, beam.createElectron()),
      mMonte(&beam, &mChamber, nullptr),
      pl(Math2::MINUS_Z_AXIS, Math2::ORIGIN_3D),
      //pls{ &pl },
      //blk(pls, 1),
      mCap0Outer(end0, radius),
      mCap1Outer(end1, radius),
      mCylOuter(end0, end1, radius),
      //shapes{ &mCylOuter, &mCap0Outer, &mCap1Outer },
      //mPillOuter(shapes, 3),
      bmm1(mat1),
      bmm2(mat2),
      r1(&mChamber, &bmm1, &blk),
      r2(&r1, &bmm2, &mPillOuter)
   {
      blk.addPlane(&pl);
      mPillOuter.addShape(&mCylOuter);
      mPillOuter.addShape(&mCap0Outer);
      mPillOuter.addShape(&mCap1Outer);
   }

   static void assertTrue(bool expr)
   {
      if (!expr) printf("false!\n");
   }

   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);;
   }

   void SumShapeTest::pointInside(double res[])
   {
      switch (Random::randomInt(3)) {
      case 0: {
         const double r = radius * Random::random();
         const double th = Random::random() * Math2::PI;
         const double phi = 2.0 * Random::random() * Math2::PI;
         res[0] = -length + r * ::cos(phi) * ::sin(th);
         res[1] = r * ::sin(phi) * ::sin(th);
         res[2] = 2.0 * radius + r * ::cos(th);
         assertTrue(mCap0Outer.contains(res));
      }
              break;
      case 1: {
         const double phi = 2.0 * Random::random() * Math2::PI;
         const double r = radius * Random::random();
         res[0] = 2.0 * length * (Random::random() - 0.5);
         res[1] = r * ::cos(phi);
         res[2] = 2.0 * radius + ::sin(phi) * r;
         assertTrue(mCylOuter.contains(res));
      }
              break;
      case 2: {
         const double r = radius * Random::random();
         const double th = Random::random() * Math2::PI;
         const double phi = 2.0 * Random::random() * Math2::PI;
         res[0] = length + r * ::cos(phi) * ::sin(th);
         res[1] = r * ::sin(phi) * ::sin(th);
         res[2] = 2.0 * radius + r * ::cos(th);
         assertTrue(mCap1Outer.contains(res));
      }
              break;
      }
   }

   static void pointOutside(double res[])
   {
      const double phi = 2.0 * Random::random() * Math2::PI;
      const double r = radius * (1.00001 + 0.9 * Random::random());
      res[0] = 3.0 * length * (Random::random() - 0.5);
      res[1] = r * ::cos(phi);
      res[2] = 2.0 * radius + ::sin(phi) * r;
   }

   void SumShapeTest::testGetFirstIntersection()
   {
      for (int i = 0; i < 1000; ++i) {
         double inside1[3];
         pointInside(inside1);
         assertTrue(mPillOuter.contains(inside1));
         double inside2[3];
         pointInside(inside2);
         assertTrue(mPillOuter.contains(inside2));
         double outside[3];
         pointOutside(outside);
         assertTrue(!mPillOuter.contains(outside));
         {
            const double t = mPillOuter.getFirstIntersection(inside1, inside2);
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
            assertTrue(mPillOuter.contains(inside1));
            assertTrue(!mPillOuter.contains(outside));
            const double t = mPillOuter.getFirstIntersection(inside1, outside);
            //if (vrml != null) {
            //   if (wr != null)
            //      wr.println("# inside to outside");
            //   vrml.drawLine(inside1, Math2::pointBetween(inside1, outside, t));
            //}
            double ptbtw0[3];
            Math2::pointBetween3d(inside1, outside, 0.99 * t, ptbtw0);
            assertTrue(mPillOuter.contains(ptbtw0));
            double ptbtw1[3];
            Math2::pointBetween3d(inside1, outside, 1.01 * t, ptbtw1);
            assertTrue(!mPillOuter.contains(ptbtw1));
            const double tp = mPillOuter.getFirstIntersection(outside, inside1);
            //if (vrml != null) {
            //   wr.println("# outside to inside");
            //   vrml.drawLine(outside, Math2::pointBetween(outside, inside1, tp));
            //   wr.flush();
            //}
            double ptbtw2[3];
            Math2::pointBetween3d(outside, inside1, 0.99 * tp, ptbtw2);
            assertTrue(!mPillOuter.contains(ptbtw2));
            double ptbtw3[3];
            Math2::pointBetween3d(outside, inside1, 1.01 * tp, ptbtw3);
            assertTrue(mPillOuter.contains(ptbtw3));
            assertTrue(t < 1.0);
            assertTrue(tp < 1.0);
            assertEquals(1.0, t + tp, 1.0e-6);
         }
         {
            const double t = mPillOuter.getFirstIntersection(inside2, outside);
            const double tp = mPillOuter.getFirstIntersection(outside, inside2);
            assertTrue(t < 1.0);
            assertTrue(tp < 1.0);
            assertEquals(1.0, t + tp, 1.0e-6);
         }
      }

      printf("SumShapeTest::testGetFirstIntersection() completed.\n");
   }

   void SumShapeTest::testAll()
   {
      mMonte.runMultipleTrajectories(1000);

      printf("SumShapeTest::testAll() completed.\n");
   }
}
