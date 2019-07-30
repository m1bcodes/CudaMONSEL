#ifndef _FROM_SI_CUH_
#define _FROM_SI_CUH_

#include <cuda_runtime.h>

namespace FromSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ extern const float KEV;
   __constant__ extern const float EV;
   __constant__ extern const float GRAM;
   __constant__ extern const float CM;

   __constant__ extern const float MICROMETER;
   //__constant__ extern const double AMU;

   __constant__ extern const float ANGSTROM;
   __constant__ extern const float TORR;

   __constant__ extern const float NANO;
   __constant__ extern const float PICO;
#else
   extern const float KEV;
   extern const float EV;
   extern const float GRAM;
   extern const float CM;

   extern const float MICROMETER;
   //extern const double AMU;

   extern const float ANGSTROM;
   extern const float TORR;

   extern const float NANO;
   extern const float PICO;
#endif

   double Torr(double pascal);
   __host__ __device__ float keV(const float e);
   __host__ __device__ float eV(const float e);
   double AMU(double kg);
   double dyne(double f);
   double gPerCC(double d);
   double angstrom(double a);
   double sqrAngstrom(double a2);
   double cmSqrPerg(double x);
   double cm(double x);
   double micrometer(double x);
   double nanometer(double x);
   double centigrade(double c);
   double fahrenheit(double k);
}

#endif