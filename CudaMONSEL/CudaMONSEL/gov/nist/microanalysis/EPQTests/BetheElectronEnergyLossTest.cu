#include "gov\nist\microanalysis\EPQTests\BetheElectronEnergyLossTest.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BetheElectronEnergyLoss.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

namespace BetheElectronEnergyLossTest
{
   static void assertEquals(double a, double b, double err)
   {
      if (::abs(a - b) >= err) printf("not equal\n");
   }

   void testOne()
   {
      assertEquals(BetheElectronEnergyLoss::JoyLuo1989.compute(Element::Mn, ToSI::keV(10.0)), -2.15748e-13, 1e-16);
      assertEquals(BetheElectronEnergyLoss::JoyLuo1989.compute(Element::U, ToSI::keV(1.0)), -3.994744E-13, 1e-16);
   }
}