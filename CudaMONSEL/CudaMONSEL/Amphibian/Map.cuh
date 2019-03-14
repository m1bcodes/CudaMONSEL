#ifndef _MAP_CUH_
#define _MAP_CUH_

#include "Set.cuh"
#include "LinkedList.cuh"
#include "Hasher.cuh"

#include <cuda_runtime.h>

#include <stdio.h>

namespace Map
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int NUM_BUCKETS = 23;
#else
   const int NUM_BUCKETS = 23;
#endif

   template<typename K, typename V>
   class Iterator;

   template<typename K, typename V>
   class Map
   {
   friend class Iterator<K, V>;
   public:
      typedef bool(*pCmp)(K, K);
      __host__ __device__ Map(Hasher::pHasher, pCmp);
      __host__ __device__ Map(const Map<K, V>&);
      __host__ __device__ Map<K, V>& operator=(const Map<K, V>&);
      __host__ __device__ ~Map();

      __host__ __device__ void Initialize();
      __host__ __device__ void ClearAndCopy(const Map<K, V>&);
      __host__ __device__ void Put(K, V);
      __host__ __device__ bool ContainsKey(K);
      __host__ __device__ V GetValue(K);
      __host__ __device__ Set::Set<K> GetKeys();
      __host__ __device__ unsigned int Hash(K);
      __host__ __device__ unsigned int HashCode();
      __host__ __device__ void DeepCopy(const Map& other);
      __host__ __device__ V Remove(K k);
      __host__ __device__ void RemoveAll();
      __host__ __device__ bool IsEmpty();
      __host__ __device__ int Size();
      __host__ __device__ V Aggregate(V(*fcn)(V));
      //__host__ __device__ LinkedListKV::Node<K, V>* AsList();
      __host__ __device__ LinkedListKV::Node<K, V>* GetBucket(int); // DEBUGGING PURPOSES

   private:
      __host__ __device__ int unsigned GetBucketIdx(K k);

      LinkedListKV::Node<K, V>* buckets[NUM_BUCKETS];
      Hasher::pHasher hasher = NULL;
      pCmp cmp = NULL;
   };

   template<typename K, typename V>
   __host__ __device__ Map<K, V>::Map(Hasher::pHasher hasher, Map<K, V>::pCmp cmp) : hasher(hasher), cmp(cmp)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename K, typename V>
   __host__ __device__ Map<K, V>::Map(const Map<K, V>& m)
   {
      //printf("called cc\n");
      ClearAndCopy(m);
   }

   template<typename K, typename V>
   __host__ __device__ Map<K, V>& Map<K, V>::operator=(const Map<K, V>& m)
   {
      //printf("called =\n");
      ClearAndCopy(m);
      return *this;
   }

   template<typename K, typename V>
   __host__ __device__ Map<K, V>::~Map()
   {
      RemoveAll();
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::Initialize()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
      hasher = NULL;
      cmp = NULL;
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::DeepCopy(const Map<K,V>& other)
   {
      RemoveAll();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         auto itr = other.buckets[k];
         while (itr != NULL) {
            Put(itr->GetKey(), itr->GetValue());
            itr = itr->GetNext();
         }
      }
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::ClearAndCopy(const Map<K, V>& other)
   {
      Initialize();
      hasher = other.hasher;
      cmp = other.cmp;

      DeepCopy(other);
   }

   template<typename K, typename V>
   __host__ __device__ V Map<K, V>::GetValue(K k)
   {
      return LinkedListKV::GetValue<K, V>(buckets[GetBucketIdx(k)], k, cmp);
   }

   template<typename K, typename V>
   __host__ __device__ Set::Set<K> Map<K, V>::GetKeys()
   {
      Set::Set<K> res(hasher, cmp);
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         auto itr = buckets[k];
         while (itr != NULL) {
            res.Put(itr->GetKey());
            itr = itr->GetNext();
         }
      }
      return res;
   }

   template<typename K, typename V>
   __host__ __device__ V Map<K, V>::Remove(K k)
   {
      return LinkedListKV::Remove<K, V>(&buckets[GetBucketIdx(k)], k, cmp);
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::RemoveAll()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKV::RemoveAll<K, V>(&buckets[k]);
      }
   }

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::Put(K k, V v)
   {
      if (!ContainsKey(k)) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(k)], k, v);
      }
      else {
         LinkedListKV::Node<K, V>* bucketItr = buckets[GetBucketIdx(k)];
         while (bucketItr != NULL) {
            if (cmp(bucketItr->GetKey(), k)) {
               bucketItr->MapVal(v, [](V a, V b) { return a; });
               break;
            }
            bucketItr = bucketItr->GetNext();
         }
      }
   }

   template<typename K, typename V>
   __host__ __device__ bool Map<K, V>::ContainsKey(K k)
   {
      return LinkedListKV::ContainsKey<K, V>(buckets[GetBucketIdx(k)], k, cmp);
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

   template<typename K, typename V>
   __host__ __device__ unsigned int Map<K, V>::GetBucketIdx(K v)
   {
      return Hash(v) % NUM_BUCKETS;
   }

   template<typename K, typename V>
   __host__ __device__ unsigned int Map<K, V>::HashCode()
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedListKV::HashCode(buckets[k], hasher);
      }

      return res;
   }

   //template<typename K, typename V>
   //__host__ __device__ LinkedListKV::Node<K, V>* Map<K, V>::AsList()
   //{
   //   LinkedListKV::Node<K, V>* res = NULL;
   //   for (int k = 0; k < NUM_BUCKETS; ++k) {
   //      auto itr = buckets[k];
   //      while (itr != NULL) {
   //         LinkedListKV::InsertHead<K, V>(&res, itr->GetKey(), itr->GetValue());
   //         itr = itr->GetNext();
   //      }
   //   }
   //   return res;
   //}

   template<typename K, typename V>
   __host__ __device__ bool Map<K, V>::IsEmpty()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != NULL) {
            return false;
         }
      }
      return true;
   }

   template<typename K, typename V>
   __host__ __device__ int Map<K, V>::Size()
   {
      int c = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         auto itr = buckets[k];
         while (itr != NULL) {
            ++c;
            itr = itr->GetNext();
         }
      }
      return c;
   }

   template<typename K, typename V>
   __host__ __device__ LinkedListKV::Node<K, V>* Map<K, V>::GetBucket(int n)
   {
      return buckets[n];
   }

   template<typename K, typename V>
   __host__ __device__ V Map<K, V>::Aggregate(V(*fcn)(V))
   {
      V res = NULL;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         auto itr = buckets[k];
         while (itr != NULL) {
            auto v = itr->GetValue();
            res += fcn(v);
            itr = itr->GetNext();
         }
      }
      return res;
   }

   template<typename K, typename V>
   class Iterator
   {
   public:
      __host__ __device__ Iterator(Map<K, V>&);
      __host__ __device__ void Reset();
      __host__ __device__ void Next();

      __host__ __device__ bool HasNext();

      __host__ __device__ K GetKey();
      __host__ __device__ V GetValue();

   private:
      LinkedListKV::Node<K, V>* ptr;
      Map<K, V>& refMap;
      int bucket;
   };

   template<typename K, typename V>
   __host__ __device__ Iterator<K, V>::Iterator(Map<K, V>& m) : refMap(m), ptr(NULL), bucket(-1)
   {
      Reset();
   }

   template<typename K, typename V>
   __host__ __device__ void Iterator<K, V>::Reset()
   {
      if (refMap.IsEmpty()) {
         ptr = NULL;
         bucket = -1;
         return;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (refMap.buckets[k] != NULL) {
            bucket = k;
            ptr = refMap.buckets[bucket];
            break;
         }
      }
   }

   template<typename K, typename V>
   __host__ __device__ void Iterator<K, V>::Next()
   {
      if (bucket == -1) {
         return;
      }
      if (ptr != NULL) {
         ptr = ptr->GetNext();
      }
      if (ptr == NULL) {
         for (int k = bucket + 1; k < NUM_BUCKETS; ++k) {
            if (refMap.buckets[k] != NULL) {
               bucket = k;
               ptr = refMap.buckets[bucket];
               return;
            }
         }
         bucket = -1;
      }
   }

   template<typename K, typename V>
   __host__ __device__ bool Iterator<K, V>::HasNext()
   {
      return bucket != -1;
   }

   template<typename K, typename V>
   __host__ __device__ K Iterator<K, V>::GetKey()
   {
      return (bucket == -1) ? NULL : ptr->GetKey();
   }

   template<typename K, typename V>
   __host__ __device__ V Iterator<K, V>::GetValue()
   {
      return (bucket == -1) ? NULL : ptr->GetValue();
   }
}

#endif
