//#ifndef _MAP_TEST_CUH_
//#define _MAP_TEST_CUH_
//
//#include <cuda_runtime.h>
//
//#include "Amphibian/Hasher.cuh"
//#include "Amphibian/Map.cuh"
//
//namespace MapTest
//{
//   class MapTest
//   {
//   public:
//      __device__ MapTest();
//      __device__ void TestInteger();
//      __device__ void TestString();
//      __device__ void TestMapOfMap();
//      __device__ void TestAggregate();
//
//      template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//      __device__ Map::Map<K, V, KCompare, VCompare, KHasher, VHasher> CreateMapA(K& k, V& v)
//      {
//         Map::Map<K, V, KCompare, VCompare, KHasher, VHasher> m1;
//         m1.Put(k, v);
//
//         return m1;
//      }
//   };
//}
//
//#endif
