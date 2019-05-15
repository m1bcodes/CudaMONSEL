#ifndef _STRATEGY_CUH_
#define _STRATEGY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include <map>
#include <set>

namespace Strategy
{
   typedef std::map<StringT, const AlgorithmClassT*> AlgorithmMap;
   typedef std::set<StringT> AlgorithmNames;
   typedef std::set<const AlgorithmClassT*> Algorithms;

   class Strategy
   {
   public:
      const AlgorithmMap& getStrategyMap() const;
      const AlgorithmNames listAlgorithmClasses() const;

      Strategy();

      Strategy& operator=(const Strategy& strat);

      void apply(const Strategy& st);
      void addAll(const Strategy& st);

      void addAlgorithm(StringT cls, const AlgorithmClassT* value);
      const AlgorithmClassT* getAlgorithm(StringT cls) const;

      Algorithms getAlgorithms() const;

      void clear();
      bool empty() const;

   private:
      AlgorithmMap mMap;
   };
}

#endif