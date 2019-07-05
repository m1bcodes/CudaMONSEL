#include "ToSI.cuh"
#include "PhysicalConstants.cuh"

namespace ToSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const double MEV = 1.60217653e-19 * 1.0e6; // PhysicalConstants::ElectronCharge * 1.0e6;
   __constant__ const double KEV = 1.60217653e-19 * 1.0e3; // PhysicalConstants::ElectronCharge * 1.0e3;
   __constant__ const double EV = 1.60217653e-19; // PhysicalConstants::ElectronCharge;
   __constant__ const double GRAM = 1.0e-3;
   __constant__ const double CM = 1.0e-2;
   __constant__ const double MICROMETER = 1.0e-6;
   // const double AMU = PhysicalConstants::UnifiedAtomicMass;
   __constant__ const double ANGSTROM = 1.0e-10;
   __constant__ const double BARN = 1.0e-28;
   __constant__ const double TORR = 101325.0 / 760; // PhysicalConstants::StandardAtmosphere / 760.0;
   __constant__ const double NANO = 1.0e-9;
   __constant__ const double PICO = 1.0e-12;
#else
   const double MEV = 1.60217653e-19 * 1.0e6; // PhysicalConstants::ElectronCharge * 1.0e6;
   const double KEV = 1.60217653e-19 * 1.0e3; // PhysicalConstants::ElectronCharge * 1.0e3;
   const double EV = 1.60217653e-19; // PhysicalConstants::ElectronCharge;
   const double GRAM = 1.0e-3;
   const double CM = 1.0e-2;
   const double MICROMETER = 1.0e-6;
   // const double AMU = PhysicalConstants::UnifiedAtomicMass;
   const double ANGSTROM = 1.0e-10;
   const double BARN = 1.0e-28;
   const double TORR = 101325.0 / 760; // PhysicalConstants::StandardAtmosphere / 760.0;
   const double NANO = 1.0e-9;
   const double PICO = 1.0e-12;
#endif

   double Torr(double torr)
   {
      return TORR * torr;
   }

   double MeV(double e)
   {
      return MEV * e;
   }

   __host__ __device__ double keV(double e)
   {
      return KEV * e;
   }

   __host__ __device__ double eV(double e)
   {
      return EV * e;
   }

   double AMU(double amu)
   {
      return PhysicalConstants::UnifiedAtomicMass * amu;
   }

   double dyne(double f)
   {
      return 1.0e-5 * f;
   }

   double gPerCC(double d)
   {
      return d * GRAM / (CM * CM * CM);
   }

   double inverse_gPerCC(double d)
   {
      return d * (CM * CM * CM) / GRAM;
   }

   double percm(double d)
   {
      return d * (1 / CM);
   }

   double angstrom(double a)
   {
      return ANGSTROM * a;
   }

   double micrometer(double m)
   {
      return MICROMETER * m;
   }

   double sqrAngstrom(double a2)
   {
      return ANGSTROM * ANGSTROM * a2;
   }

   double barn(double a2)
   {
      return BARN * a2; // meters^2 per barn
   }

   double cmSqrPerg(double x)
   {
      return ((CM * CM) / GRAM) * x;
   }

   double cm(double x)
   {
      return x * CM;
   }

   double g(double x)
   {
      return x * GRAM;
   }

   double centigrade(double k)
   {
      return k + PhysicalConstants::IcePoint;
   }

   double fahrenheit(double f)
   {
      return centigrade((f - 32.0) * 5.0 / 9.0);
   }
}
