#include "gov\nist\microanalysis\EPQLibrary\Strategy.cuh"

namespace Strategy
{
   __host__ __device__ const AlgorithmMap& Strategy::getStrategyMap() const
   {
      return mMap;
   }

   //const AlgorithmNames Strategy::listAlgorithmClasses() const
   //{
   //   AlgorithmNames ret;
   //   for (auto p : mMap) {
   //      ret.insert(p.first);
   //   }

   //   return ret;
   //}

   __host__ __device__ Strategy::Strategy()
   {
   }

   __host__ __device__ Strategy& Strategy::operator=(const Strategy& strat)
   {
      mMap = strat.mMap;
      return *this;
   }

   __host__ __device__ void Strategy::apply(const Strategy& st)
   {
      for (auto me : st.mMap) {
         if (mMap.find(me.first) != mMap.end()) {
            mMap.Put(me.first, me.second); // mMap[k] = me.second;
         }
      }
   }

   __host__ __device__ void Strategy::addAll(const Strategy& st)
   {
      for (auto me : st.mMap) {
         mMap.Put(me.first, me.second); // mMap[me.first] = me.second;
      }
   }

   __host__ __device__ void Strategy::addAlgorithm(StringT cls, AlgorithmClassT const * value)
   {
      //if (!cls.isAssignableFrom(value.getClass()))
      //   throw new IllegalArgumentException(value.toString() + " is not derived from " + cls.toString());
      mMap.Put(cls, value); //mMap[cls] = value;
   }

   __host__ __device__ const AlgorithmClassT* Strategy::getAlgorithm(StringT cls) const
   {
      //return mMap.at(cls);
      auto itr = mMap.find(cls);
      if (itr != mMap.end()) return itr->second;
      else return nullptr;
   }

   __host__ __device__ Algorithms Strategy::getAlgorithms() const
   {
      Algorithms ret;
      for (auto e : mMap) {
         ret.insert(e.second);
      }
      return ret;
   }

   __host__ __device__ void Strategy::clear()
   {
      mMap.clear();
   }

   __host__ __device__ bool Strategy::empty() const
   {
      return mMap.empty();
   }
}