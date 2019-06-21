// file: gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh

#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"

#include <cuda_runtime.h>

namespace ZeroCSD
{
   class ZeroCSD : public SlowingDownAlgT
   {
   public:
      __host__ __device__ ZeroCSD();

      double compute(double d, const ElectronT* pe) const override;
      void setMaterial(const SEmaterialT* mat) override;

      StringT toString() const override;
   };
}