#ifndef _ToSI_CUH_
#define _ToSI_CUH_

#include <cuda_runtime.h>

namespace ToSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ extern const float MEV;
   __constant__ extern const float KEV;
   __constant__ extern const float EV;
   __constant__ extern const float GRAM;
   __constant__ extern const float CM;
   __constant__ extern const float MICROMETER;
   //__constant__extern const double AMU;
   __constant__ extern const float ANGSTROM;
   __constant__ extern const float BARN;
   __constant__ extern const float TORR;
   __constant__ extern const float GIGA;
   __constant__ extern const float NANO;
   __constant__ extern const float PICO;
#else
   extern const float MEV;
   extern const float KEV;
   extern const float EV;
   extern const float GRAM;
   extern const float CM;
   extern const float MICROMETER;
   //extern const double AMU;
   extern const float ANGSTROM;
   extern const float BARN;
   extern const float TORR;
   extern const float GIGA;
   extern const float NANO;
   extern const float PICO;
#endif

   double Torr(double torr);
   double MeV(double e);
   __host__ __device__ float keV(float e);
   __host__ __device__ float eV(float e);
   __host__ __device__ float AMU(float amu);
   double dyne(double f);
   double gPerCC(double d);
   double inverse_gPerCC(double d);
   double percm(double d);
   __host__ __device__ float angstrom(float a);
   double micrometer(double m);
   __host__ __device__ float sqrAngstrom(float a2);
   double barn(double a2);
   double cmSqrPerg(double x);
   double cm(double x);
   double g(double x);
   double centigrade(double k);
   double fahrenheit(double f);
};

#endif
