#include "gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh"

namespace ZeroCSD
{
   __host__ __device__ ZeroCSD::ZeroCSD()
   {
   }

   __host__ __device__ double ZeroCSD::compute(double d, const ElectronT* pe) const
   {
      return 0.;
   }

   __host__ __device__ void ZeroCSD::setMaterial(const SEmaterialT* mat)
   {
   }

   __host__ __device__ StringT ZeroCSD::toString() const
   {
      return "ZeroCSD()";
   }
}