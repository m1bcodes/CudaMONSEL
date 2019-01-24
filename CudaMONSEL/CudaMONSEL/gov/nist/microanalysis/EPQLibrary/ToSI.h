#ifndef ToSI_H
#define ToSI_H

#include "PhysicalConstants.h"

class ToSI
{
public:
   static const double MEV;
   static const double KEV;
   static const double EV;
   static const double GRAM;
   static const double CM;
   static const double MICROMETER;
   //static const double AMU;
   static const double ANGSTROM;
   static const double BARN;
   static const double TORR;

   static const double NANO;
   static const double PICO;

   static double Torr(double torr)
   {
      return TORR * torr;
   }

   static double MeV(double e)
   {
      return MEV * e;
   }

   static double keV(double e)
   {
      return KEV * e;
   }

   static double eV(double e)
   {
      return EV * e;
   }

   static double AMU(double amu)
   {
      return PhysicalConstants::UnifiedAtomicMass * amu;
   }

   static double dyne(double f)
   {
      return 1.0e-5 * f;
   }

   static double gPerCC(double d)
   {
      return d * GRAM / (CM * CM * CM);
   }

   static double inverse_gPerCC(double d)
   {
      return d * (CM * CM * CM) / GRAM;
   }

   static double percm(double d)
   {
      return d * (1 / CM);
   }

   static double angstrom(double a)
   {
      return ANGSTROM * a;
   }

   static double micrometer(double m)
   {
      return MICROMETER * m;
   }

   static double sqrAngstrom(double a2)
   {
      return ANGSTROM * ANGSTROM * a2;
   }

   static double barn(double a2)
   {
      return BARN * a2; // meters^2 per barn
   }

   static double cmSqrPerg(double x)
   {
      return ((CM * CM) / GRAM) * x;
   }

   static double cm(double x)
   {
      return x * CM;
   }

   static double g(double x)
   {
      return x * GRAM;
   }

   static double centigrade(double k)
   {
      return k + PhysicalConstants::IcePoint;
   }

   static double fahrenheit(double f)
   {
      return centigrade((f - 32.0) * 5.0 / 9.0);
   }
};

#endif
