#ifndef _MAP_CUH_
#define _MAP_CUH_

#include "Set.cuh"
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
   class MapIterator;

   template<typename K, typename V>
   class Map
   {
   friend class MapIterator<K, V>;
   public:
      typedef bool(*pCmp)(K, K);
      __host__ __device__ Map(Hasher::pHasher, pCmp);
      //__host__ __device__ Map(const Map&);
      //__host__ __device__ ~Map();
      __host__ __device__ void Put(K, V);
      __host__ __device__ bool ContainsKey(K);
      __host__ __device__ V GetValue(K);
      __host__ __device__ Set::Set<K> GetKeys();
      __host__ __device__ unsigned int Hash(K);
      __host__ __device__ unsigned int HashCode();
      __host__ __device__ void DeepCopy(const Map& other);
      __host__ __device__ V Remove(K k);
      __host__ __device__ void RemoveAll();
      __host__ __device__ LinkedListKV::Node<K, V>* AsList();
      __host__ __device__ bool IsEmpty();
      __host__ __device__ int Size();
      __host__ __device__ V Aggregate(V (*fcn)(V));
      __host__ __device__ LinkedListKV::Node<K, V>* GetBucket(int); // DEBUGGING PURPOSES
      //__host__ __device__ LinkedListKV::Node<K, V>* GetBucket(int n);

   private:
      __host__ __device__ int unsigned GetBucketIdx(K k);

      LinkedListKV::Node<K, V>* buckets[NUM_BUCKETS];
      Hasher::pHasher hasher;
      pCmp cmp;
   };

   template<typename K, typename V>
   __host__ __device__ Map<K, V>::Map(Hasher::pHasher hasher, Map<K, V>::pCmp cmp) : hasher(hasher), cmp(cmp)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   //template<typename K, typename V>
   //__host__ __device__ Map<K, V>::Map(const Map<K, V>& m)
   //{
   //   DeepCopy(m);
   //}

   //template<typename K, typename V>
   //__host__ __device__ Map<K, V>::~Map()
   //{
   //   RemoveAll();
   //}

   template<typename K, typename V>
   __host__ __device__ void Map<K, V>::DeepCopy(const Map& other)
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

   template<typename K, typename V>
   __host__ __device__ LinkedListKV::Node<K, V>* Map<K, V>::AsList()
   {
      LinkedListKV::Node<K, V>* res = NULL;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         auto itr = buckets[k];
         while (itr != NULL) {
            LinkedListKV::InsertHead<K, V>(&res, itr->GetKey(), itr->GetValue());
            itr = itr->GetNext();
         }
      }
      return res;
   }

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
   class MapIterator
   {
   public:
      __host__ __device__ MapIterator(const Map<K, V>&);
      __host__ __device__ void Reset();
      __host__ __device__ void Next();

      __host__ __device__ K GetKey();
      __host__ __device__ V GetValue();

   private:
      LinkedListKV::Node<K, V>* ptr;
      Map<K, V> refMap;
      int bucket;
   };

   template<typename K, typename V>
   __host__ __device__ MapIterator<K, V>::MapIterator(const Map<K, V>& m) : refMap(m), ptr(NULL), bucket(0)
   {
   }

   template<typename K, typename V>
   __host__ __device__ void MapIterator<K, V>::Reset()
   {
      bucket = 0;
      ptr = refMap.buckets[bucket];
   }

   template<typename K, typename V>
   __host__ __device__ void MapIterator<K, V>::Next()
   {
      if (ptr != NULL) {
         ptr = ptr->GetNext();
      }
      if (ptr == NULL) {
         bucket = (bucket + 1) % NUM_BUCKETS;
         ptr = refMap.buckets[bucket];
      }
   }

   template<typename K, typename V>
   __host__ __device__ K MapIterator<K, V>::GetKey()
   {
      return ptr->GetKey();
   }

   template<typename K, typename V>
   __host__ __device__ V MapIterator<K, V>::GetValue()
   {
      return ptr->GetValue();
   }
}

#endif
