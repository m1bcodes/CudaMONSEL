#include "FromSI.cuh"
#include "PhysicalConstants.cuh"

namespace FromSI
{
   __device__ const double KEV = 1.0 / (1.60217653e-19 * 1.0e3); // 1.0 / (PhysicalConstants::ElectronCharge * 1.0e3);
   __device__ const double EV = 1.0 / 1.60217653e-19; // 1.0 / PhysicalConstants::ElectronCharge;
   __device__ const double GRAM = 1000.0;
   __device__ const double CM = 100.0;

   __device__ const double MICROMETER = 1.0e6;
   //__device__ const double AMU = 1.0 / PhysicalConstants.UnifiedAtomicMass;

   __device__ const double ANGSTROM = 1.0e10;
   __device__ const double TORR = 760.0 / 101325.0; // 760.0 / PhysicalConstants::StandardAtmosphere;

   __device__ const double NANO = 1.0e9;
   __device__ const double PICO = 1.0e12;

   __device__ double Torr(double pascal) {
      return TORR * pascal;
   }

   __device__ double keV(double e) {
      return e * KEV;
   }

   __device__ double eV(double e) {
      return e * EV;
   }

   __device__ double AMU(double kg) {
      return 1.0 / PhysicalConstants::UnifiedAtomicMass * kg;
   }

   __device__ double dyne(double f) {
      return 1.0e5 * f;
   }

   __device__ double gPerCC(double d) {
      return d * (GRAM / (CM * CM * CM));
   }

   __device__ double angstrom(double a) {
      return ANGSTROM * a;
   }

   __device__ double sqrAngstrom(double a2) {
      return ANGSTROM * ANGSTROM * a2;
   }

   __device__ double cmSqrPerg(double x) {
      return ((CM * CM) / GRAM) * x;
   }

   __device__ double cm(double x) {
      return x * CM;
   }

   __device__ double micrometer(double x) {
      return x * MICROMETER;
   }

   __device__ double nanometer(double x) {
      return x * NANO;
   }

   __device__ double centigrade(double c) {
      return c - PhysicalConstants::IcePoint;
   }

   __device__ double fahrenheit(double k) {
      return 32.0 + centigrade(k) * 9.0 / 5.0;
   }
}