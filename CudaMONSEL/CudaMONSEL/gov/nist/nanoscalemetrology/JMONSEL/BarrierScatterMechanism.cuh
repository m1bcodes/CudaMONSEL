// file: gov\nist\nanoscalemetrology\JMONSEL\BarrierScatterMechanism.cuh

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace BarrierScatterMechanism
{
   class BarrierScatterMechanism
   {
   public:
      __host__ __device__ virtual ElectronT* barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const = 0;
   };
}
