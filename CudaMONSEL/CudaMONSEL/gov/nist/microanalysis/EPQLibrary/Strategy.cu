#include "gov\nist\microanalysis\EPQLibrary\Strategy.cuh"

namespace Strategy
{
   const AlgorithmMap& Strategy::getStrategyMap() const
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

   Strategy::Strategy()
   {
   }

   Strategy& Strategy::operator=(const Strategy& strat)
   {
      mMap = strat.mMap;
      return *this;
   }

   void Strategy::apply(const Strategy& st)
   {
      for (auto me : st.mMap) {
         auto k = me.first;
         if (mMap.find(k) != mMap.end()) {
            mMap[k] = me.second;
         }
      }
   }

   void Strategy::addAll(const Strategy& st)
   {
      for (auto me : st.mMap) {
         mMap[me.first] = me.second;
      }
   }

   void Strategy::addAlgorithm(StringT cls, AlgorithmClassT const * const value)
   {
      //if (!cls.isAssignableFrom(value.getClass()))
      //   throw new IllegalArgumentException(value.toString() + " is not derived from " + cls.toString());
      mMap[cls] = value;
   }

   const AlgorithmClassT* Strategy::getAlgorithm(StringT cls) const
   {
      //return mMap.at(cls);
      auto itr = mMap.find(cls);
      return (itr != mMap.end()) ? itr->second : nullptr;
   }

   Algorithms Strategy::getAlgorithms() const
   {
      Algorithms ret;
      for (auto e : mMap) {
         ret.insert(e.second);
      }
      return ret;
   }

   void Strategy::clear()
   {
      mMap.clear();
   }

   bool Strategy::empty() const
   {
      return mMap.empty();
   }
}