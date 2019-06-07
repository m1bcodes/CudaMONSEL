#ifndef _UNORDERED_MAP_CUH_
#define _UNORDERED_MAP_CUH_

#include "Amphibian/unordered_set.cuh"
#include "Amphibian/LinkedList.cuh"
#include "Amphibian/Hasher.cuh"

#include <cuda_runtime.h>

#include <stdio.h>

namespace Map
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static const int NUM_BUCKETS = 23;
#else
   static const int NUM_BUCKETS = 23;
#endif

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   class Map
   {
      typedef LinkedListKV::Node<K, V>* LinkedListKVPtr;

   public:
      __host__ __device__ Map();
      __host__ __device__ Map(const Map&);
      __host__ __device__ Map(int);
      __host__ __device__ Map& operator=(const Map&);
      __host__ __device__ ~Map();

      __host__ __device__ bool operator==(const Map&) const;

      __host__ __device__ void Initialize();
      __host__ __device__ void ClearAndCopy(const Map&);
      __host__ __device__ void Put(const K&, const V&);
      __host__ __device__ bool ContainsKey(const K&) const;
      __host__ __device__ LinkedListKVPtr Find(const K& k) const;
      __host__ __device__ bool GetValue(const K&, V&) const;
      __host__ __device__ amp::unordered_set<K, KCompare, KHasher> GetKeys();
      __host__ __device__ unsigned int Hash(const K&) const;
      __host__ __device__ unsigned int HashCode() const;
      __host__ __device__ void DeepCopy(const Map& other);
      __host__ __device__ bool Remove(const K&, V&);
      __host__ __device__ bool Remove(const K&);
      __host__ __device__ void RemoveAll();
      __host__ __device__ bool IsEmpty() const;
      __host__ __device__ int Size() const;
      __host__ __device__ double Aggregate(double(*fcn)(double)) const;
      //__host__ __device__ LinkedListKV::Node<K, V>* AsList();

      //__host__ __device__ Hasher::pHasher GetHasher(); // DEBUGGING PURPOSES
      //__host__ __device__ pKeyCmp GetKeyCmp(); // DEBUGGING PURPOSES
      //__host__ __device__ pValCmp GetValCmp(); // DEBUGGING PURPOSES

      class Iterator
      {
      public:
         __host__ __device__ Iterator(Map<K, V, KCompare, VCompare, KHasher, VHasher>&);
         __host__ __device__ void Reset();
         __host__ __device__ void Next();

         __host__ __device__ bool HasNext();

         __host__ __device__ K& GetKey() const;
         __host__ __device__ V& GetValue() const;

      private:
         LinkedListKV::Node<K, V>* ptr;
         Map<K, V, KCompare, VCompare, KHasher, VHasher>& refMap;
         int bucket;
      };

   private:
      __host__ __device__ int unsigned GetBucketIdx(const K& k) const;
      __host__ __device__ LinkedListKVPtr GetBucket(int); // DEBUGGING PURPOSES

      LinkedListKVPtr buckets[NUM_BUCKETS];
      KCompare kcmp;
      VCompare vcmp;
      KHasher khasher;
      VHasher vhasher;
   };

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ Hasher::pHasher Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetHasher()
   //{
   //   return hasher;
   //}

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::pKeyCmp Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeyCmp()
   //{
   //   return kcmp;
   //}

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::pValCmp Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValCmp()
   //{
   //   return vcmp;
   //}

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
   {
      //printf("called cc\n");
      ClearAndCopy(m);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map(int a)
   {
      Initialize();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>& Map<K, V, KCompare, VCompare, KHasher, VHasher>::operator=(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
   {
      //printf("called =\n");
      if (&m == this) {
         return *this;
      }
      ClearAndCopy(m);
      return *this;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::~Map()
   {
      //printf("called des\n");
      RemoveAll();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::operator==(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& other) const
   {
      if (this == &other) return true;

      if (Size() != other.Size()) {
         return false;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = other.buckets[k];
         while (itr != NULL) {
            K k1 = itr->GetKey();
            V v0, v1;
            bool found0 = GetValue(k1, v0);
            bool found1 = other.GetValue(k1, v1);
            if (found0 && !found1 || !found0 && found1) { // only one of the keys is NULL
               return false;
            }
            if (found0 && found1) { // both values are not NULL
               if (!vcmp(v0, v1)) { // values are different
                  return false;
               }
            }
            itr = itr->GetNext();
         }
      }
      return true;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Initialize()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::DeepCopy(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
   {
      RemoveAll();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         LinkedListKVPtr itr = other.buckets[k];
         while (itr != NULL) {
            K k = itr->GetKey();
            V v = itr->GetValue();
            Put(k, v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::ClearAndCopy(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
   {
      Initialize();

      DeepCopy(other);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValue(const K& k, V& v) const
   {
      //LinkedListKV::GetValue<K, V>(buckets[0], k, kcmp);
      return LinkedListKV::GetValue<K, V, KCompare>(buckets[GetBucketIdx(k)], k, v);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ amp::unordered_set<K, KCompare, KHasher> Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeys()
   {
      amp::unordered_set<K, KCompare, KHasher> res;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != NULL) {
            res.Put(itr->GetKey());
            itr = itr->GetNext();
         }
      }
      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::Remove(const K& k, V& ret)
   {
      return LinkedListKV::Remove<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, ret);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::Remove(const K& k)
   {
      V v;
      return Remove(k, v);
      //return LinkedListKV::Remove<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, v);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::RemoveAll()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKV::RemoveAll<K, V>(&buckets[k]);
      }
   }

   template<typename V>
   __host__ __device__ static const V& GetFirstParam(const V& a, V& b)
   {
      return a;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Put(const K& k, const V& v)
   {
      LinkedListKVPtr itr = Find(k);
      if (itr == nullptr) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(k)], k, v);
      }
      else {
         //LinkedListKVPtr bucketItr = buckets[GetBucketIdx(k)];
         //while (bucketItr != NULL) {
         //   const K& tmpK = bucketItr->GetKey();
         //   if (kcmp(tmpK, k)) {
         //      bucketItr->MapVal(v, GetFirstParam<V>);
         //      break;
         //   }
         //   bucketItr = bucketItr->GetNext();
         //}
         itr->MapVal(v, GetFirstParam<V>);
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::ContainsKey(const K& k) const
   {
      return LinkedListKV::ContainsKey<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::LinkedListKVPtr Map<K, V, KCompare, VCompare, KHasher, VHasher>::Find(const K& k) const
   {
      return LinkedListKV::Find<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::Hash(const K& k) const
   {
      return khasher(k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucketIdx(const K& v) const
   {
      return Hash(v) % NUM_BUCKETS;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::HashCode() const
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedListKV::HashCode<K, V, KHasher, VHasher>(buckets[k]);
      }

      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::IsEmpty() const
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != NULL) {
            return false;
         }
      }
      return true;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ int Map<K, V, KCompare, VCompare, KHasher, VHasher>::Size() const
   {
      int c = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != NULL) {
            ++c;
            itr = itr->GetNext();
         }
      }
      return c;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::LinkedListKVPtr Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucket(int n)
   {
      return buckets[n];
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ double Map<K, V, KCompare, VCompare, KHasher, VHasher>::Aggregate(double(*fcn)(double)) const
   {
      double res = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != NULL) {
            const V& v = itr->GetValue();
            res += fcn(v);
            itr = itr->GetNext();
         }
      }
      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool AreEqual(Map<K, V, KCompare, VCompare, KHasher, VHasher>& a, Map<K, V, KCompare, VCompare, KHasher, VHasher>& b)
   {
      if (&a == &b) return true;
      return a == b;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::Iterator(Map<K, V, KCompare, VCompare, KHasher, VHasher>& m) : refMap(m), ptr(NULL), bucket(-1)
   {
      Reset();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::Reset()
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

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::Next()
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

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::HasNext()
   {
      return bucket != -1;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ K& Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::GetKey() const
   {
      if (bucket == -1 || ptr == NULL) {
         printf("Illegal call to map iterator GetKey(): NULL pointer, no more element");
      }
      return ptr->GetKey();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ V& Map<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator::GetValue() const
   {
      if (bucket == -1 || ptr == NULL) {
         printf("Illegal call to map iterator GetValue(): NULL pointer, no more element");
      }
      return ptr->GetValue();
   }
}

#endif
