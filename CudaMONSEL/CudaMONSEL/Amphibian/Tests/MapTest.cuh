//#ifndef _MAP_TEST_CUH_
//#define _MAP_TEST_CUH_
//
//#include <cuda_runtime.h>
//
//#include "..\Hasher.cuh"
//#include "..\Map.cuh"
//
//namespace MapTest
//{
//   class MapTest
//   {
//   public:
//      __device__ MapTest();
//      __device__ void TestInteger();
//      __device__ void TestString();
//
//      template<typename K, typename V>
//      __device__ Map::Map<K, V> CreateMapA(bool(*kcmp)(K, K), bool(*vcmp)(V, V), K k, V v)
//      {
//         Map::Map<K, V> m1(DefaultHasher, kcmp, vcmp);
//         m1.Put(k, v);
//
//         return m1;
//      }
//
//   private:
//      Hasher::pHasher DefaultHasher;
//   };
//}
//
//#endif