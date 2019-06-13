#include "gov\nist\microanalysis\EPQLibrary\Gas.cuh"
#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"

namespace Gas
{
   static long serialVersionUID = 0x43;

   static const double RYDBERG_CONSTANT = PhysicalConstants::BoltzmannConstant * PhysicalConstants::AvagadroNumber; // J
   // per
   // K
   // per
   // mole
   static const double MOLAR_VOLUME = 0.022414; // cubic meters per
   // mole at 0 C & 1 ATM
   static const double STANDARD_PRESSURE = 101325.0; // pascals
   static const double STANDARD_TEMPERATURE = 273.15; // kelvin

   double Gas::getMassPerSubunit()
   {
      double res = 0.0;
      for (auto elm : getElementSet())
         res += atomicPercent(*elm) * elm->getMass();
      return res;
   }

   Gas::Gas() : Material(0.0), mTemperature(STANDARD_TEMPERATURE), mPressure(STANDARD_PRESSURE)
   {
   }

   Gas::Gas(ElementT elms[], int stoic[], int len, double pressure, double temperature, const char* name) : Material(0.0), mTemperature(temperature), mPressure(pressure)
   {
      for (int i = 0; i < len; ++i)
         addElementByStoiciometry(elms[i], stoic[i]);
      
      // n = (P*V)/(R*T)
      double n = (1.0 /* m^3 */* pressure) / (RYDBERG_CONSTANT * temperature);
      setDensity((getMassPerSubunit() * PhysicalConstants::AvagadroNumber) * n);
      setName(name);
   }

   void Gas::setTemperature(double newTemp)
   {
      if (newTemp != mTemperature) {
         mPressure *= (newTemp / mTemperature);
         mTemperature = newTemp;
      }
   }

   double Gas::getTemperature()
   {
      return mTemperature;
   }

   double Gas::getPressure()
   {
      return mPressure;
   }

   Gas::Gas(const Gas& other) : MaterialT(other), mPressure(other.mPressure), mTemperature(other.mTemperature)
   {
   }

   void Gas::replicate(const Gas& gas)
   {
      MaterialT::replicate(gas);
      mPressure = gas.mPressure;
      mTemperature = gas.mTemperature;
   }

   Gas Gas::clone()
   {
      Gas res;
      replicate(res);
      return res;
   }

   //private void writeObject(java.io.ObjectOutputStream out)
   //   throws IOException{
   //   out.writeDouble(mPressure);
   //   out.writeDouble(mTemperature);
   //}

   //   private void readObject(java.io.ObjectInputStream in)
   //   throws IOException, ClassNotFoundException{
   //   mPressure = in.readDouble();
   //   mTemperature = in.readDouble();
   //}
}
