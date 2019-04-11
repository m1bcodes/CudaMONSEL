#ifndef _FROM_SI_CUH_
#define _FROM_SI_CUH_

#include <cuda_runtime.h>

namespace FromSI
{
   extern const double KEV;
   extern const double EV;
   extern const double GRAM;
   extern const double CM;

   extern const double MICROMETER;
   //extern const double AMU;

   extern const double ANGSTROM;
   extern const double TORR;

   extern const double NANO;
   extern const double PICO;

   double Torr(double pascal);
   double keV(double e);
   double eV(double e);
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