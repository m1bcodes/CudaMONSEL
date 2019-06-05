#ifndef _GANACHAUD_MOKRANI_PHONON_INELASTIC_SM_CUH_
#define _GANACHAUD_MOKRANI_PHONON_INELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace GanachaudMokraniPhononInelasticSM
{
   class GanachaudMokraniPhononInelasticSM : public ScatterMechanismT
   {
   public:
      GanachaudMokraniPhononInelasticSM(double ratemultiplier, double phononE, double temperature, double eps0, double epsInfinity);

      double scatterRate(const ElectronT* pe) override;
      ElectronT* scatter(ElectronT* pe) override;
      void setMaterial(const MaterialT* mat) override;

      const char* toString();

   private:
      const double ratemultiplier;
      const double phononE; // Energy of the phonon mode
      const double occupationFactor; // Typically ~ 1/2 at room temperature
      const double epsRatio; // (eps0-epsinfinity)/esp0/epsinfinity
      const double prefactor;
      const double temperature;
      const double eps0;
      const double epsInfinity;

      char buff[100]; // work around for toString
   };
}

#endif