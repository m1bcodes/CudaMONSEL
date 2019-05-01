#ifndef _RADOMIZED_SCATTER_FACTORY_CUH_
#define _RADOMIZED_SCATTER_FACTORY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmClass.cuh"

namespace RandomizedScatterFactory
{
   class RandomizedScatterFactory : public AlgorithmClassT
   {
   public:
      AlgorithmClassT const * const * getAllImplementations() const override;

      virtual const RandomizedScatterT& get(const ElementT& elm) const = 0;

   protected:
      RandomizedScatterFactory(StringT name, const ReferenceT& ref);
   };
}

#endif