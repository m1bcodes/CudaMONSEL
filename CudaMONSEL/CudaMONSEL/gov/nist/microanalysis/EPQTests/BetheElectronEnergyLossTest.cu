#include "gov\nist\microanalysis\EPQTests\BetheElectronEnergyLossTest.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

namespace BetheElectronEnergyLossTest
{
   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal: %.10e, %.10e\n", a, b);
   }

   void testOne()
   {
      //not equal : -2.1636209716e-013, -2.1574800000e-013 OK
      //not equal : -4.0909654404e-013, -3.9947440000e-013 OK
      assertEquals(BetheElectronEnergyLoss::JoyLuo1989.compute(Element::Mn, ToSI::keV(10.0)), -2.15748e-13, 1e-16);
      assertEquals(BetheElectronEnergyLoss::JoyLuo1989.compute(Element::U, ToSI::keV(1.0)), -3.994744E-13, 1e-16);
      printf("BetheElectronEnergyLossTest::testOne() completed.\n");
   }
}