#ifndef _ALGORITHM_CLASS_CUH_
#define _ALGORITHM_CLASS_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\AlgorithmUser.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Reference.cuh"

#include <cuda_runtime.h>

namespace AlgorithmClass
{
   class AlgorithmClass : public AlgorithmUserT
   {
   protected:
      __host__ __device__ AlgorithmClass(StringT, StringT, const ReferenceT&);
      //AlgorithmClass(StringT clss, StringT name, StringT ref);

   public:
      virtual AlgorithmClassT const * const * getAllImplementations() const = 0;

      int compareTo(const AlgorithmClass&) const;
      __host__ __device__ const StringT toString() const;
      const StringT getAlgorithmClass() const;
      __host__ __device__ const StringT getName() const;
      const StringT getReference() const;
      __host__ __device__ const ReferenceT& getReferenceObj() const;

   private:
      const StringT mClass;
      const StringT mName;
      const ReferenceT& mReference;
   };
};

#endif