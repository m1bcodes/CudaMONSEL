#ifndef _ToSI_CUH_
#define _ToSI_CUH_

#include <cuda_runtime.h>

namespace ToSI
{
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

   double Torr(double torr);
   double MeV(double e);
   double keV(double e);
   double eV(double e);
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
