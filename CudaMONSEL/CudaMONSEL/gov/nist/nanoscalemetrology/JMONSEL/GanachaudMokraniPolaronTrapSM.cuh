#ifndef _GANACHAUD_MOKIRANI_POLARON_TRAP_SM_CUH_
#define _GANACHAUD_MOKIRANI_POLARON_TRAP_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace GanachaudMokraniPolaronTrapSM
{
   class GanachaudMokraniPolaronTrapSM : public ScatterMechanismT
   {
   public:
      __host__ __device__ GanachaudMokraniPolaronTrapSM(double prefactor, double extinction);

      __host__ __device__ double scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      double prefactor;
      double extinction;
      double CUTOFF;
   };
}

#endif