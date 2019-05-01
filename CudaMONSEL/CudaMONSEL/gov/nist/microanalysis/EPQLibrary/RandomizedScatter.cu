#include "gov\nist\microanalysis\EPQLibrary\RandomizedScatter.cuh"

namespace RandomizedScatter
{
   RandomizedScatter::RandomizedScatter(StringT name, const ReferenceT& ref) : AlgorithmClassT("Elastic cross-section", name, ref)
   {
   }
   
   void RandomizedScatter::initializeDefaultStrategy()
   {
   }

   AlgorithmClassT const * const * RandomizedScatter::getAllImplementations() const
   {
      return NULL;
   }
}
