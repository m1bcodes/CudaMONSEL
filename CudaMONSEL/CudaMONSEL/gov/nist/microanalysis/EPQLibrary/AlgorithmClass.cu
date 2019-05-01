#include "AlgorithmClass.cuh"

namespace AlgorithmClass
{
   AlgorithmClass::AlgorithmClass(StringT cls, StringT name, const ReferenceT& ref) : mClass(cls), mName(name), mReference(ref)
   {
   }

   AlgorithmClass::AlgorithmClass(StringT clss, StringT name, StringT ref) : mReference(Reference::CrudeReference(ref)), mClass(clss), mName(name)
   {
   }

   int AlgorithmClass::compareTo(const AlgorithmClass& o) const
   {
      return toString() != o.toString();
   }

   const StringT AlgorithmClass::toString() const
   {
      return mClass + "[" + mName + "]";
   }

   const StringT AlgorithmClass::getAlgorithmClass() const
   {
      return mClass;
   }

   const StringT AlgorithmClass::getName() const
   {
      return mName;
   }

   const StringT AlgorithmClass::getReference() const
   {
      return mReference.getShortForm();
   }

   const ReferenceT& AlgorithmClass::getReferenceObj() const
   {
      return mReference;
   }
}
