#include "ToSI.cuh"
#include "PhysicalConstants.cuh"

namespace ToSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const float MEV = 1.60217653e-19 * 1.0e6; // PhysicalConstants::ElectronCharge * 1.0e6;
   __constant__ const float KEV = 1.60217653e-19 * 1.0e3; // PhysicalConstants::ElectronCharge * 1.0e3;
   __constant__ const float EV = 1.60217653e-19; // PhysicalConstants::ElectronCharge;
   __constant__ const float GRAM = 1.0e-3;
   __constant__ const float CM = 1.0e-2;
   __constant__ const float MICROMETER = 1.0e-6;
   // const double AMU = PhysicalConstants::UnifiedAtomicMass;
   __constant__ const float ANGSTROM = 1.0e-10;
   __constant__ const float BARN = 1.0e-28;
   __constant__ const float TORR = 101325.0 / 760; // PhysicalConstants::StandardAtmosphere / 760.0;
   __constant__ const float GIGA = 1.0e9;
   __constant__ const float NANO = 1.0e-9;
   __constant__ const float PICO = 1.0e-12;
#else
   const float MEV = 1.60217653e-19 * 1.0e6; // PhysicalConstants::ElectronCharge * 1.0e6;
   const float KEV = 1.60217653e-19 * 1.0e3; // PhysicalConstants::ElectronCharge * 1.0e3;
   const float EV = 1.60217653e-19; // PhysicalConstants::ElectronCharge;
   const float GRAM = 1.0e-3;
   const float CM = 1.0e-2;
   const float MICROMETER = 1.0e-6;
   // const double AMU = PhysicalConstants::UnifiedAtomicMass;
   const float ANGSTROM = 1.0e-10;
   const float BARN = 1.0e-28;
   const float TORR = 101325.0 / 760; // PhysicalConstants::StandardAtmosphere / 760.0;
   const float GIGA = 1.0e9;
   const float NANO = 1.0e-9;
   const float PICO = 1.0e-12;
#endif

   double Torr(double torr)
   {
      return TORR * torr;
   }

   double MeV(double e)
   {
      return MEV * e;
   }

   __host__ __device__ float keV(float e)
   {
      return KEV * e;
   }

   __host__ __device__ float eV(float e)
   {
      return EV * e;
   }

   __host__ __device__ float AMU(float amu)
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

   __host__ __device__ float angstrom(float a)
   {
      return ANGSTROM * a;
   }

   double micrometer(double m)
   {
      return MICROMETER * m;
   }

   __host__ __device__ float sqrAngstrom(float a2)
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
