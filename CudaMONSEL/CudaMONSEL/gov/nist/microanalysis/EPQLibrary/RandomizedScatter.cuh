#ifndef _RANDOMIZED_SCATTER_CUH_
#define _RANDOMIZED_SCATTER_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace RandomizedScatter
{
   class RandomizedScatter : public AlgorithmClassT
   {
   protected:
      __host__ __device__ RandomizedScatter(StringT name, const ReferenceT& ref);

      __host__ __device__ void initializeDefaultStrategy() override;
      AlgorithmClass const * const * getAllImplementations() const override;

   public:
      __host__ __device__ virtual double randomScatteringAngle(const double energy) const = 0;
      __host__ __device__ virtual double totalCrossSection(const double energy) const = 0;
      __host__ __device__ virtual const ElementT& getElement() const = 0;
   };
}

#endif