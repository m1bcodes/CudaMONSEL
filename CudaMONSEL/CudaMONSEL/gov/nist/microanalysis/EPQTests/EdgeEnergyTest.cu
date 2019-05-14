#include "gov\nist\microanalysis\EPQTests\EdgeEnergyTest.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AtomicShell.cuh"
#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace EdgeEnergyTest
{
   static void assertEquals(double a, double b, double err)
   {
      if (abs(a - b) >= err) printf("wrong: %.10e, %.10e\n", a, b);
   }

   void testOne()
   {
      AtomicShellT sh(Element::W, AtomicShell::LIII);
      //assertEquals(EdgeEnergy::Chantler2005.compute(sh), EdgeEnergy::DTSA.compute(sh), ToSI::eV(10.0));
      //assertEquals(EdgeEnergy::NISTxrtdb.compute(sh), EdgeEnergy::DTSA.compute(sh), ToSI::eV(10.0));
      //assertEquals(EdgeEnergy::Wernish84.compute(sh), EdgeEnergy::DTSA.compute(sh), ToSI::eV(20.0));

      //assertEquals(EdgeEnergy::Chantler2005.compute(sh), EdgeEnergy::DTSA.compute(sh), ToSI::eV(10.0));
      assertEquals(EdgeEnergy::NISTxrtdb.compute(sh), EdgeEnergy::Chantler2005.compute(sh), ToSI::eV(10.0));
      assertEquals(EdgeEnergy::Wernish84.compute(sh), EdgeEnergy::Chantler2005.compute(sh), ToSI::eV(20.0));

      printf("EdgeEnergyTest::testOne() completed.\n");
   }
}
