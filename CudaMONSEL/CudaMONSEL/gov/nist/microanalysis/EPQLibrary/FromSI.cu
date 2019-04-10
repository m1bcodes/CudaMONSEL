#include "FromSI.cuh"
#include "PhysicalConstants.cuh"

namespace FromSI
{
   const double KEV = 1.0 / (1.60217653e-19 * 1.0e3); // 1.0 / (PhysicalConstants::ElectronCharge * 1.0e3);
   const double EV = 1.0 / 1.60217653e-19; // 1.0 / PhysicalConstants::ElectronCharge;
   const double GRAM = 1000.0;
   const double CM = 100.0;

   const double MICROMETER = 1.0e6;
   //const double AMU = 1.0 / PhysicalConstants.UnifiedAtomicMass;

   const double ANGSTROM = 1.0e10;
   const double TORR = 760.0 / 101325.0; // 760.0 / PhysicalConstants::StandardAtmosphere;

   const double NANO = 1.0e9;
   const double PICO = 1.0e12;

   double Torr(double pascal) {
      return TORR * pascal;
   }

   double keV(double e) {
      return e * KEV;
   }

   double eV(double e) {
      return e * EV;
   }

   double AMU(double kg) {
      return 1.0 / PhysicalConstants::UnifiedAtomicMass * kg;
   }

   double dyne(double f) {
      return 1.0e5 * f;
   }

   double gPerCC(double d) {
      return d * (GRAM / (CM * CM * CM));
   }

   double angstrom(double a) {
      return ANGSTROM * a;
   }

   double sqrAngstrom(double a2) {
      return ANGSTROM * ANGSTROM * a2;
   }

   double cmSqrPerg(double x) {
      return ((CM * CM) / GRAM) * x;
   }

   double cm(double x) {
      return x * CM;
   }

   double micrometer(double x) {
      return x * MICROMETER;
   }

   double nanometer(double x) {
      return x * NANO;
   }

   double centigrade(double c) {
      return c - PhysicalConstants::IcePoint;
   }

   double fahrenheit(double k) {
      return 32.0 + centigrade(k) * 9.0 / 5.0;
   }
}