#include "gov\nist\microanalysis\EPQTests\MonteCarloSSTest.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\SimpleBlock.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"

namespace MonteCarloSSTest
{
   static void assertEquals(int a, int b)
   {
      if (a != b) {
         printf("not equal: %d, %d\n", a, b);
      }
   }

   static void assertTrue(bool expr)
   {
      if (!expr) printf("false!\n");
   }

   void testOne()
   {
      assertEquals(MonteCarloSS::ScatterEvent, 1);

      const double center[] = {
         0.0,
         0.0,
         0.0
      };
      const double beamEnergy = ToSI::keV(20.0);
      const double beamDia = 1.0e-8;
      const GaussianBeamT beam(beamDia, beamEnergy, center);
      const double ChamberRadius = 0.1;
      SphereT sphere(center, ChamberRadius);
      NullMaterialScatterModelT NULL_MSM;
      RegionT mChamber(nullptr, &NULL_MSM, &sphere);
      MonteCarloSST mcss(&beam, &mChamber, beam.createElectron());
      // / TODO Restore these tests....
      // assertEquals(mcss.getMassAbsorptionCoefficient().compute(Element.Fe,
      // ToSI.eV(4511)), ToSI.inverse_gPerCC(ToSI.percm(187.7)),
      // ToSI.inverse_gPerCC(ToSI.percm(187.7)) / 10);
      // assertEquals(mcss.getMassAbsorptionCoefficient().compute(Element.U,
      // ToSI.eV(2697)), ToSI.inverse_gPerCC(ToSI.percm(1068.)),
      // ToSI.inverse_gPerCC(ToSI.percm(1068.)) / 5);
      // assertEquals(mcss.getMassAbsorptionCoefficient().compute(Element.Fe,
      // ToSI.eV(4511.0)), ToSI.inverse_gPerCC(ToSI.percm(187.7)),
      // ToSI.inverse_gPerCC(ToSI.percm(187.7)) / 10);
      // assertEquals(mcss.getMassAbsorptionCoefficient().compute(Element.Ti,
      // ToSI.keV(20.0),false), 1.54e-21, 0.1e-21);
      {
         const double o[] = {
            0.001,
            -0.001,
            0.0004
         };
         const double t[] = {
            1.0,
            1.0,
            1.0
         };
         double u = mcss.getChamber()->getShape2()->getFirstIntersection(o, t);
         if (!(u > 0.0)) printf("false: %.10e\n", u);
         const double res[] = {
            o[0] + u * (t[0] - o[0]),
            o[1] + u * (t[1] - o[1]),
            o[2] + u * (t[2] - o[2])
         };
         bool b = mcss.getChamber()->getShape()->contains(res);
         assertTrue(b);
      }
      {
         const double c0[] = {
            -1.0,
            -1.0,
            -1.0
         };
         const double c1[] = {
            1.0,
            1.0,
            1.0
         };
         SimpleBlockT sb(c0, c1);
         {
            const double pos0[] = {
               -1.0,
               -1.0,
               -3.0
            };
            const double pos1[] = {
               1.0,
               1.0,
               2.0
            };
            const double u = sb.getFirstIntersection(pos0, pos1);
            const double res[] = {
               pos0[0] + (pos1[0] - pos0[0]) * u,
               pos0[1] + (pos1[1] - pos0[1]) * u,
               pos0[2] + (pos1[2] - pos0[2]) * u
            };
            assertTrue(sb.contains(res));
         }
      }

      printf("MonteCarloSSTest::testOne() completed.\n");
   }
}
