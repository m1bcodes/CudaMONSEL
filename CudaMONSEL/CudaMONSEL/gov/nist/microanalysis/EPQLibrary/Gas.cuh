#ifndef _GAS_CUH_
#define _GAS_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"

namespace Gas
{
   class Gas : public MaterialT
   {
   protected:
      Gas();
      Gas(const Gas&);
      void replicate(const Gas& gas);

   public:
      double getMassPerSubunit();
      Gas(ElementT elms[], int stoic[], int len, double pressure, double temperature, const char* name);
      void setTemperature(double newTemp);
      double getTemperature();
      double getPressure();
      Gas clone();

   private:
      double mPressure; // Pressure in pascals
      double mTemperature; // Temperature in Kelvin
   };
}

#endif