// file: gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh

#ifndef _ZERO_CSD_CUH_
#define _ZERO_CSD_CUH_

#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"

#include <cuda_runtime.h>

namespace ZeroCSD
{
   class ZeroCSD : public SlowingDownAlgT
   {
   public:
      __host__ __device__ ZeroCSD();

      __host__ __device__::SlowingDownAlg::data_type compute(::SlowingDownAlg::data_type d, const ElectronT* pe) const override;
      __host__ __device__ void setMaterial(const SEmaterialT* mat) override;

      __host__ __device__ StringT toString() const override;
   };
}

#endif