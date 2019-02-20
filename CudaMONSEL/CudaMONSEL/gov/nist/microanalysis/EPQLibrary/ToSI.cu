#include "ToSI.cuh"
#include "PhysicalConstants.cuh"

namespace ToSI
{
   __device__ const double MEV = 1.60217653e-19 * 1.0e6; // PhysicalConstants::ElectronCharge * 1.0e6;
   __device__ const double KEV = 1.60217653e-19 * 1.0e3; // PhysicalConstants::ElectronCharge * 1.0e3;
   __device__ const double EV = 1.60217653e-19; // PhysicalConstants::ElectronCharge;
   __device__ const double GRAM = 1.0e-3;
   __device__ const double CM = 1.0e-2;
   __device__ const double MICROMETER = 1.0e-6;
   // __device__ const double AMU = PhysicalConstants::UnifiedAtomicMass;
   __device__ const double ANGSTROM = 1.0e-10;
   __device__ const double BARN = 1.0e-28;
   __device__ const double TORR = 101325.0 / 760; // PhysicalConstants::StandardAtmosphere / 760.0;
   __device__ const double NANO = 1.0e-9;
   __device__ const double PICO = 1.0e-12;

   __device__ double Torr(double torr)
   {
      return TORR * torr;
   }

   __device__ double MeV(double e)
   {
      return MEV * e;
   }

   __device__ double keV(double e)
   {
      return KEV * e;
   }

   __device__ double eV(double e)
   {
      return EV * e;
   }

   __device__ double AMU(double amu)
   {
      return PhysicalConstants::UnifiedAtomicMass * amu;
   }

   __device__ double dyne(double f)
   {
      return 1.0e-5 * f;
   }

   __device__ double gPerCC(double d)
   {
      return d * GRAM / (CM * CM * CM);
   }

   __device__ double inverse_gPerCC(double d)
   {
      return d * (CM * CM * CM) / GRAM;
   }

   __device__ double percm(double d)
   {
      return d * (1 / CM);
   }

   __device__ double angstrom(double a)
   {
      return ANGSTROM * a;
   }

   __device__ double micrometer(double m)
   {
      return MICROMETER * m;
   }

   __device__ double sqrAngstrom(double a2)
   {
      return ANGSTROM * ANGSTROM * a2;
   }

   __device__ double barn(double a2)
   {
      return BARN * a2; // meters^2 per barn
   }

   __device__ double cmSqrPerg(double x)
   {
      return ((CM * CM) / GRAM) * x;
   }

   __device__ double cm(double x)
   {
      return x * CM;
   }

   __device__ double g(double x)
   {
      return x * GRAM;
   }

   __device__ double centigrade(double k)
   {
      return k + PhysicalConstants::IcePoint;
   }

   __device__ double fahrenheit(double f)
   {
      return centigrade((f - 32.0) * 5.0 / 9.0);
   }
}
