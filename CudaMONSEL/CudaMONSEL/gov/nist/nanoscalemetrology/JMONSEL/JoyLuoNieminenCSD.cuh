#ifndef _JOY_LUO_NIEMINEN_CSD_CUH_
#define _JOY_LUO_NIEMINEN_CSD_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"

namespace JoyLuoNieminenCSD
{
   class JoyLuoNieminenCSD : public SlowingDownAlgT
   {
   public:
      __host__ __device__ JoyLuoNieminenCSD(SEmaterialT& mat, float breakE);
      JoyLuoNieminenCSD(MaterialT& mat, float bh, float breakE);

      void setMaterial(const MaterialT& mat, float bh);

      __host__ __device__ void init();

      __host__ __device__ void setMaterial(const SEmaterialT* mat) override;
      __host__ __device__::SlowingDownAlg::data_type compute(::SlowingDownAlg::data_type d, const ElectronT* pe) const override;

      __host__ __device__ float compute(const float len, const float kE) const;

      float getBreakE() const;
      void setBreakE(float breakE);

      __host__ __device__ StringT toString() const override;

   private:
      MaterialT& mat;

      float bh; // barrier height

      int nce; // # constituent elements

      /* Combinations of variables that we need */
      VectorXf recipJ; // 1.166/J where J=ionization energy of const.
      // elem.

      VectorXf coef; // 78500 rho c[i]Z[i]/A[i] = leading coefficient in

      VectorXf beta; // 1=1.166(wf+1eV)/J[i]

      float bhplus1eV; // BH + 1eV

      float minEforTracking;

      float breakE; // Dividing line between Joy/Luo and Nieminen
      float gamma; // The proportionality constant in Nieminen's formula
   };
}

#endif