#ifndef _MAP_TEST_CUH_
#define _MAP_TEST_CUH_

#include <cuda_runtime.h>

#include "Amphibian/Hasher.cuh"
#include "Amphibian/unordered_map.cuh"

namespace MapTest
{
   class MapTest
   {
   public:
      __host__ __device__ MapTest();
      __host__ __device__ void testInteger();
      __host__ __device__ void testString();
      __host__ __device__ void testMapOfMap();
      __host__ __device__ void testAggregate();

      template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
      __host__ __device__ amp::unordered_map<K, V, KCompare, VCompare, KHasher, VHasher> CreateMapA(K& k, V& v)
      {
         amp::unordered_map<K, V, KCompare, VCompare, KHasher, VHasher> m1;
         m1.Put(k, v);

         return m1;
      }
   };
}

#endif
