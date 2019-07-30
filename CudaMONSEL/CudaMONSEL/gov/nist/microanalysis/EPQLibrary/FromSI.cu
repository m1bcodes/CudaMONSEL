#include "FromSI.cuh"
#include "PhysicalConstants.cuh"

namespace FromSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const float KEV = 1.0 / (1.60217653e-19 * 1.0e3); // 1.0 / (PhysicalConstants::ElectronCharge * 1.0e3);
   __constant__ const float EV = 1.0 / 1.60217653e-19; // 1.0 / PhysicalConstants::ElectronCharge;
   __constant__ const float GRAM = 1000.0;
   __constant__ const float CM = 100.0;

   __constant__ const float MICROMETER = 1.0e6;
   //const double AMU = 1.0 / PhysicalConstants.UnifiedAtomicMass;

   __constant__ const float ANGSTROM = 1.0e10;
   __constant__ const float TORR = 760.0 / 101325.0; // 760.0 / PhysicalConstants::StandardAtmosphere;

   __constant__ const float NANO = 1.0e9;
   __constant__ const float PICO = 1.0e12;
#else
   const float KEV = 1.0 / (1.60217653e-19 * 1.0e3); // 1.0 / (PhysicalConstants::ElectronCharge * 1.0e3);
   const float EV = 1.0 / 1.60217653e-19; // 1.0 / PhysicalConstants::ElectronCharge;
   const float GRAM = 1000.0;
   const float CM = 100.0;

   const float MICROMETER = 1.0e6;
   //const double AMU = 1.0 / PhysicalConstants.UnifiedAtomicMass;

   const float ANGSTROM = 1.0e10;
   const float TORR = 760.0 / 101325.0; // 760.0 / PhysicalConstants::StandardAtmosphere;

   const float NANO = 1.0e9;
   const float PICO = 1.0e12;
#endif

   double Torr(double pascal) {
      return TORR * pascal;
   }

   __host__ __device__ float keV(const float e) {
      return e * KEV;
   }

   __host__ __device__ float eV(const float e) {
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