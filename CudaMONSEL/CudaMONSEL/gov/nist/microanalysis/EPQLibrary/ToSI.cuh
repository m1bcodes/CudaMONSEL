#ifndef _ToSI_CUH_
#define _ToSI_CUH_

#include <cuda_runtime.h>

namespace ToSI
{
   extern __device__ const double MEV;
   extern __device__ const double KEV;
   extern __device__ const double EV;
   extern __device__ const double GRAM;
   extern __device__ const double CM;
   extern __device__ const double MICROMETER;
   //extern __device__ const double AMU;
   extern __device__ const double ANGSTROM;
   extern __device__ const double BARN;
   extern __device__ const double TORR;
   extern __device__ const double NANO;
   extern __device__ const double PICO;

   __device__ double Torr(double torr);
   __device__ double MeV(double e);
   __device__ double keV(double e);
   __device__ double eV(double e);
   __device__ double AMU(double amu);
   __device__ double dyne(double f);
   __device__ double gPerCC(double d);
   __device__ double inverse_gPerCC(double d);
   __device__ double percm(double d);
   __device__ double angstrom(double a);
   __device__ double micrometer(double m);
   __device__ double sqrAngstrom(double a2);
   __device__ double barn(double a2);
   __device__ double cmSqrPerg(double x);
   __device__ double cm(double x);
   __device__ double g(double x);
   __device__ double centigrade(double k);
   __device__ double fahrenheit(double f);
};

#endif
