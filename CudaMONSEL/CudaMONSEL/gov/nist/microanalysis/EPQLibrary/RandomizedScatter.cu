#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"

namespace RandomizedScatter
{
   __host__ __device__ RandomizedScatter::RandomizedScatter(StringT name, const ReferenceT& ref) : AlgorithmClassT("Elastic cross-section", name, ref)
   {
   }
   
   __host__ __device__ void RandomizedScatter::initializeDefaultStrategy()
   {
   }

   AlgorithmClassT const * const * RandomizedScatter::getAllImplementations() const
   {
      return nullptr;
   }
}
