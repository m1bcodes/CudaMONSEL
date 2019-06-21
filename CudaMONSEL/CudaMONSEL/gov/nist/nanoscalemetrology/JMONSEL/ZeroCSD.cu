#include "gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh"

namespace ZeroCSD
{
   __host__ __device__ ZeroCSD::ZeroCSD()
   {
   }

   double ZeroCSD::compute(double d, const ElectronT* pe) const
   {
      return 0.;
   }

   void ZeroCSD::setMaterial(const SEmaterialT* mat)
   {
   }

   StringT ZeroCSD::toString() const
   {
      return "ZeroCSD()";
   }
}