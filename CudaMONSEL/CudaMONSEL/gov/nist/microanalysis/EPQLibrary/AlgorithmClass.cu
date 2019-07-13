#include "AlgorithmClass.cuh"

namespace AlgorithmClass
{
   __host__ __device__ AlgorithmClass::AlgorithmClass(StringT cls, StringT name, const ReferenceT& ref) : mClass(cls), mName(name), mReference(ref)
   {
   }

   //AlgorithmClass::AlgorithmClass(StringT clss, StringT name, StringT ref) : mReference(Reference::CrudeReference(ref)), mClass(clss), mName(name)
   //{
   //}

   int AlgorithmClass::compareTo(const AlgorithmClass& o) const
   {
      return toString() != o.toString();
   }

   __host__ __device__ const StringT AlgorithmClass::toString() const
   {
      return mClass + " [" + mName + "]";
   }

   const StringT AlgorithmClass::getAlgorithmClass() const
   {
      return mClass;
   }

   __host__ __device__ const StringT AlgorithmClass::getName() const
   {
      return mName;
   }

   const StringT AlgorithmClass::getReference() const
   {
      return mReference.getShortForm();
   }

   __host__ __device__ const ReferenceT& AlgorithmClass::getReferenceObj() const
   {
      return mReference;
   }
}
