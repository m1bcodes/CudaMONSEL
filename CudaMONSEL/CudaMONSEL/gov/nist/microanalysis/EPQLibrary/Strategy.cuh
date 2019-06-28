#ifndef _STRATEGY_CUH_
#define _STRATEGY_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include <map>
#include <set>

#include "Amphibian\unordered_set.cuh"
#include "Amphibian\unordered_map.cuh"
#include "Amphibian\String.cuh"
#include "Amphibian\Hasher.cuh"

namespace Strategy
{
   //typedef std::map<StringT, const AlgorithmClassT*> AlgorithmMap;
   ////typedef std::set<StringT> AlgorithmNames;
   //typedef std::set<const AlgorithmClassT*> Algorithms;

   struct CompareFcn
   {
      __host__ __device__ inline bool operator() (const AlgorithmClassT* lhs, const AlgorithmClassT* rhs) const
      {
         return lhs == rhs;
      }
   };

   struct HashFcn
   {
      __host__ __device__ inline unsigned int operator() (const AlgorithmClassT* ptr) const
      {
         return Hasher::Hash((char *)ptr, sizeof(ptr));
      }
   };

   typedef amp::unordered_map<StringT, const AlgorithmClassT*, amp::string_cmp, CompareFcn, amp::string_hash, HashFcn> AlgorithmMap;
   //typedef amp::unordered_set<StringT, amp::string_hash, amp::string_cmp> AlgorithmNames;
   typedef amp::unordered_set<const AlgorithmClassT*, HashFcn, CompareFcn> Algorithms;

   class Strategy
   {
   public:
      __host__ __device__ const AlgorithmMap& getStrategyMap() const;
      //const AlgorithmNames listAlgorithmClasses() const;

      __host__ __device__ Strategy();

      __host__ __device__ Strategy& operator=(const Strategy& strat);

      __host__ __device__ void apply(const Strategy& st);
      __host__ __device__ void addAll(const Strategy& st);

      __host__ __device__ void addAlgorithm(StringT cls, const AlgorithmClassT* value);
      __host__ __device__ const AlgorithmClassT* getAlgorithm(StringT cls) const;

      __host__ __device__ Algorithms getAlgorithms() const;

      __host__ __device__ void clear();
      __host__ __device__ bool empty() const;

   private:
      AlgorithmMap mMap;
   };
}

#endif