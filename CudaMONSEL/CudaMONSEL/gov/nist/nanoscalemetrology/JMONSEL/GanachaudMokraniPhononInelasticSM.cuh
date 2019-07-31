#ifndef _GANACHAUD_MOKRANI_PHONON_INELASTIC_SM_CUH_
#define _GANACHAUD_MOKRANI_PHONON_INELASTIC_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace GanachaudMokraniPhononInelasticSM
{
   typedef ::ScatterMechanism::data_type data_type;

   class GanachaudMokraniPhononInelasticSM : public ScatterMechanismT
   {
   public:
      __host__ __device__ GanachaudMokraniPhononInelasticSM(data_type ratemultiplier, data_type phononE, data_type temperature, data_type eps0, data_type epsInfinity);

      __host__ __device__ data_type scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

      const char* toString();

   private:
      const data_type ratemultiplier;
      const data_type phononE; // Energy of the phonon mode
      const data_type occupationFactor; // Typically ~ 1/2 at room temperature
      const data_type epsRatio; // (eps0-epsinfinity)/esp0/epsinfinity
      const data_type prefactor;
      const data_type temperature;
      const data_type eps0;
      const data_type epsInfinity;

      char buff[100]; // work around for toString
   };
}

#endif