#include "gov\nist\microanalysis\EPQTests\AtomicShellTest.cuh"

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AtomicShell.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"

namespace AtomicShellTest
{
   static void assertTrue(bool expr)
   {
      if (!expr) printf("not true\n");
   }

   static void assertTrue(StringT a, StringT b)
   {
      if (a != b) printf("not equal: %s, %s\n", a.c_str(), b.c_str());
   }

   static void assertTrue(int a, int b)
   {
      if (a != b) printf("not equal: %d, %d\n", a, b);
   }

   static void assertTrue(double a, double b)
   {
      if (a != b) printf("not equal: %d, %d\n", a, b);
   }

   static void assertEquals(double a, double b, double err)
   {
      if (abs(a - b) > err) printf("not equal: %.10e, %.10e\n", a, b);
   }

   static void assertTrue(const ElementT& a, const ElementT& b)
   {
      if (!(a == b)) printf("not equal: %s, %s\n", a.toAbbrev(), b.toAbbrev());
   }

   void testOne()
   {
      assertTrue(AtomicShell::getAtomicName(AtomicShell::LIII), "2P3/2");
      assertTrue(AtomicShell::getFamilyName(AtomicShell::getFamily(AtomicShell::LIII)), ("L"));
      assertTrue(AtomicShell::getIUPACName(AtomicShell::LIII), ("L3"));
      assertTrue(AtomicShell::getSiegbahnName(AtomicShell::LIII), ("LIII"));
      assertTrue(AtomicShell::getPrincipalQuantumNumber(AtomicShell::LIII), 2);

      AtomicShellT as(Element::Cu, AtomicShell::LIII);
      assertTrue(as.getElement(), Element::Cu);
      assertTrue(as.getCapacity(), 4);
      assertEquals(as.getEdgeEnergy(), ToSI::keV(0.933), 0.010);
      assertTrue(as.getElement(), Element::Cu);
      assertEquals(FromSI::eV(1), FromSI::EV, 10.0);
      assertEquals(FromSI::eV(as.getEnergy()), 933.0, 10.0);
      assertTrue(as.getFamily(), AtomicShell::LFamily);
      assertTrue(as.getShell(), AtomicShell::LIII);
      assertTrue(as.getOrbitalAngularMomentum(), 1);
      assertTrue(as.getTotalAngularMomentum(), 1.5);
      assertTrue(AtomicShell::isValid(AtomicShell::LIII));
      assertTrue(as.getGroundStateOccupancy(), 4);

      printf("AtomicShell::testOne() completed.\n");
   }
}
