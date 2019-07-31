#ifndef _GANACHAUD_MOKIRANI_POLARON_TRAP_SM_CUH_
#define _GANACHAUD_MOKIRANI_POLARON_TRAP_SM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ScatterMechanism.cuh"

namespace GanachaudMokraniPolaronTrapSM
{
   typedef ::ScatterMechanism::data_type data_type;

   class GanachaudMokraniPolaronTrapSM : public ScatterMechanismT
   {
   public:
      __host__ __device__ GanachaudMokraniPolaronTrapSM(data_type prefactor, data_type extinction);

      __host__ __device__ data_type scatterRate(const ElectronT* pe) override;
      __host__ __device__ ElectronT* scatter(ElectronT* pe) override;
      __host__ __device__ void setMaterial(const MaterialT* mat) override;

      StringT toString() const;

   private:
      data_type prefactor;
      data_type extinction;
      data_type CUTOFF;
   };
}

#endif