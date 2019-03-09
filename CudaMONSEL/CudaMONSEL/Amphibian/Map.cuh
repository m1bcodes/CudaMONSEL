#ifndef _MAP_CUH_
#define _MAP_CUH_

#include "LinkedList.cuh"
#include "Hasher.cuh"

#include <cuda_runtime.h>

namespace Map
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int NUM_BUCKETS = 23;
#else
   const int NUM_BUCKETS = 23;
#endif

   template<typename K, typename V>
   class Map
   {
   public:
      typedef bool(*pCmp)(K, K);
      __host__ __device__ Map(Hasher::pHasher, pCmp);
      __host__ __device__ void Put(K, V);
      __host__ __device__ bool ContainsKey(K);
      __host__ __device__ unsigned int Hash(K);
      __host__ __device__ LinkedListKV::Node<K, V>* GetBucket(int n);

   private:
      LinkedListKV::Node<K, V>* buckets[NUM_BUCKETS] = { NULL };
      Hasher::pHasher hasher;
      pCmp cmp;
   };

   template<typename K, typename V>
   __host__ __device__ Map<K, V>::Map(Hasher::pHasher hasher, Map<K, V>::pCmp cmp) : hasher(hasher), cmp(cmp)
   {
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::Put(K k, V v)
   {
      auto bucketIdx = Hash(v) % NUM_BUCKETS;
      if (!ContainsKey(k)) {
         LinkedListKV::InsertHead<K, V>(&buckets[bucketIdx], k, v);
      }
   }

   template<typename K, typename V>
   __host__ __device__ bool Map<K, V>::ContainsKey(K k)
   {
      auto bucketIdx = Hash(v) % NUM_BUCKETS;
      return LinkedListKV::ContainsKey<K, V>(buckets[bucketIdx], k, cmp);
   }

   template<typename K, typename V>
   __host__ __device__ unsigned int Map<K, V>::Hash(K v)
   {
      return hasher((char*)&v, sizeof(v));
   }

   //template<typename T>
   //__host__ __device__ LinkedList::Node<T>* Set<T>::GetBucket(int n)
   //{
   //   return buckets[n % NUM_BUCKETS];
   //}
}

#endif
