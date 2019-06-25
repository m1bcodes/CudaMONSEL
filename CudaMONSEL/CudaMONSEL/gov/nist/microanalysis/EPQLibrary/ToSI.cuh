#ifndef _ToSI_CUH_
#define _ToSI_CUH_

#include <cuda_runtime.h>

namespace ToSI
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ extern const double MEV;
   __constant__ extern const double KEV;
   __constant__ extern const double EV;
   __constant__ extern const double GRAM;
   __constant__ extern const double CM;
   __constant__ extern const double MICROMETER;
   //__constant__extern const double AMU;
   __constant__ extern const double ANGSTROM;
   __constant__ extern const double BARN;
   __constant__ extern const double TORR;
   __constant__ extern const double NANO;
   __constant__ extern const double PICO;
#else
   extern const double MEV;
   extern const double KEV;
   extern const double EV;
   extern const double GRAM;
   extern const double CM;
   extern const double MICROMETER;
   //extern const double AMU;
   extern const double ANGSTROM;
   extern const double BARN;
   extern const double TORR;
   extern const double NANO;
   extern const double PICO;
#endif

   double Torr(double torr);
   double MeV(double e);
   double keV(double e);
   __host__ __device__ double eV(double e);
   double AMU(double amu);
   double dyne(double f);
   double gPerCC(double d);
   double inverse_gPerCC(double d);
   double percm(double d);
   double angstrom(double a);
   double micrometer(double m);
   double sqrAngstrom(double a2);
   double barn(double a2);
   double cmSqrPerg(double x);
   double cm(double x);
   double g(double x);
   double centigrade(double k);
   double fahrenheit(double f);
};

#endif
