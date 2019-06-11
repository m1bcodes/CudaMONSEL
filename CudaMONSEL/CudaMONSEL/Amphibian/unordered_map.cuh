#ifndef _UNORDERED_MAP_CUH_
#define _UNORDERED_MAP_CUH_

#include "Amphibian/unordered_set.cuh"
#include "Amphibian/LinkedList.cuh"
#include "Amphibian/Hasher.cuh"

#include <cuda_runtime.h>

#include <stdio.h>

namespace unordered_map
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static const int NUM_BUCKETS = 23;
#else
   static const int NUM_BUCKETS = 23;
#endif

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   class unordered_map
   {
      typedef LinkedListKV::Node<K, V>* LinkedListKVPtr;

   public:
      __host__ __device__ unordered_map();
      __host__ __device__ unordered_map(const unordered_map&);
      __host__ __device__ unordered_map(int);
      __host__ __device__ unordered_map& operator=(const unordered_map&);
      __host__ __device__ ~unordered_map();

      __host__ __device__ bool operator==(const unordered_map&) const;

      // capacity
      __host__ __device__ bool empty() const;
      __host__ __device__ int size() const;

      // modifiers
      __host__ __device__ void insert(LinkedListKVPtr);
      __host__ __device__ bool erase(const K&, V&);
      __host__ __device__ bool erase(const K&);

      // element lookup
      __host__ __device__ LinkedListKVPtr find(const K& k) const;

      // hash
      __host__ __device__ unsigned int hashCode() const;

      __host__ __device__ void Initialize();
      __host__ __device__ void ClearAndCopy(const unordered_map&);
      __host__ __device__ void Put(const K&, const V&);
      __host__ __device__ bool ContainsKey(const K&) const;
      __host__ __device__ bool GetValue(const K&, V&) const;
      __host__ __device__ amp::unordered_set<K, KCompare, KHasher> GetKeys();
      __host__ __device__ unsigned int Hash(const K&) const;
      __host__ __device__ void DeepCopy(const unordered_map& other);
      __host__ __device__ void RemoveAll();
      __host__ __device__ double Aggregate(double(*fcn)(double)) const;
      //__host__ __device__ LinkedListKV::Node<K, V>* AsList();

      //__host__ __device__ Hasher::pHasher GetHasher(); // DEBUGGING PURPOSES
      //__host__ __device__ pKeyCmp GetKeyCmp(); // DEBUGGING PURPOSES
      //__host__ __device__ pValCmp GetValCmp(); // DEBUGGING PURPOSES

      class iterator
      {
      public:
         __host__ __device__ iterator(unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>&);
         __host__ __device__ void Reset();
         __host__ __device__ void Next();
         __host__ __device__ void End();

         __host__ __device__ bool HasNext();

         __host__ __device__ void operator++();
         __host__ __device__ bool operator!=(const iterator&) const;
         __host__ __device__ bool operator==(const iterator&) const;
         __host__ __device__ LinkedListKV::Node<K, V>& operator*();
         __host__ __device__ void operator=(const iterator&);

         __host__ __device__ K& GetKey() const;
         __host__ __device__ V& GetValue() const;

         //class
         //{
         //public:
         //   __host__ __device__ iterator & operator = (const iterator &i) { return value = i; }
         //   __host__ __device__ operator iterator() const { return itr.GetKey(); }

         //private:
         //   iterator itr;
         //} first;

         //class
         //{
         //public:
         //   __host__ __device__ iterator & operator = (const iterator &i) { return value = i; }
         //   __host__ __device__ operator iterator() const { return itr.GetValue(); }

         //private:
         //   iterator itr;
         //} second;

      private:
         LinkedListKV::Node<K, V>* ptr;
         unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& refMap;
         int bucket;
      };

      // iterators
      __host__ __device__ iterator begin();
      __host__ __device__ iterator end();
      __host__ __device__ iterator cbegin();
      __host__ __device__ iterator cend();

   private:
      __host__ __device__ int unsigned GetBucketIdx(const K& k) const;
      __host__ __device__ LinkedListKVPtr GetBucket(int); // DEBUGGING PURPOSES

      LinkedListKVPtr buckets[NUM_BUCKETS];
   };

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ Hasher::pHasher unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetHasher()
   //{
   //   return hasher;
   //}

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::pKeyCmp unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeyCmp()
   //{
   //   return kcmp;
   //}

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::pValCmp unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValCmp()
   //{
   //   return vcmp;
   //}

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::unordered_map()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = nullptr;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::unordered_map(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
   {
      //printf("called cc\n");
      ClearAndCopy(m);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::unordered_map(int a)
   {
      Initialize();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::operator=(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
   {
      //printf("called =\n");
      if (&m != this) ClearAndCopy(m);
      return *this;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::~unordered_map()
   {
      //printf("called des\n");
      RemoveAll();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::operator==(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& other) const
   {
      VCompare vcmp;
      if (this == &other) return true;

      if (size() != other.size()) {
         return false;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = other.buckets[k];
         while (itr != nullptr) {
            K k1 = itr->GetKey();
            V v0, v1;
            bool found0 = GetValue(k1, v0);
            bool found1 = other.GetValue(k1, v1);
            if (found0 && !found1 || !found0 && found1) { // only one of the keys is nullptr
               return false;
            }
            if (found0 && found1) { // both values are not nullptr
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
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::Initialize()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = nullptr;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::DeepCopy(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
   {
      RemoveAll();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         LinkedListKVPtr itr = other.buckets[k];
         while (itr != nullptr) {
            K k = itr->GetKey();
            V v = itr->GetValue();
            Put(k, v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::ClearAndCopy(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
   {
      Initialize();

      DeepCopy(other);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValue(const K& k, V& v) const
   {
      //LinkedListKV::GetValue<K, V>(buckets[0], k, kcmp);
      return LinkedListKV::GetValue<K, V, KCompare>(buckets[GetBucketIdx(k)], k, v);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ amp::unordered_set<K, KCompare, KHasher> unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeys()
   {
      amp::unordered_set<K, KCompare, KHasher> res;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != nullptr) {
            res.Put(itr->GetKey());
            itr = itr->GetNext();
         }
      }
      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::erase(const K& k, V& ret)
   {
      return LinkedListKV::Remove<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, ret);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::erase(const K& k)
   {
      V v;
      return erase(k, v);
      //return LinkedListKV::erase<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, v);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::RemoveAll()
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
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::insert(LinkedListKV::Node<K, V>* node)
   {
      LinkedListKVPtr itr = find(node->GetKey());
      if (itr == nullptr) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(node->GetKey())], node);
      }
      else {
         itr->MapVal(node->GetValue(), GetFirstParam<V>);
         delete node;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::Put(const K& k, const V& v)
   {
      LinkedListKVPtr itr = find(k);
      if (itr == nullptr) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(k)], k, v);
      }
      else {
         //LinkedListKVPtr bucketItr = buckets[GetBucketIdx(k)];
         //while (bucketItr != nullptr) {
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
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::ContainsKey(const K& k) const
   {
      return LinkedListKV::ContainsKey<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::LinkedListKVPtr unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::find(const K& k) const
   {
      return LinkedListKV::Find<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::Hash(const K& k) const
   {
      KHasher khasher;
      return khasher(k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucketIdx(const K& v) const
   {
      return Hash(v) % NUM_BUCKETS;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::hashCode() const
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedListKV::HashCode<K, V, KHasher, VHasher>(buckets[k]);
      }

      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::empty() const
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != nullptr) {
            return false;
         }
      }
      return true;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::size() const
   {
      int c = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != nullptr) {
            ++c;
            itr = itr->GetNext();
         }
      }
      return c;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::LinkedListKVPtr unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucket(int n)
   {
      return buckets[n];
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ double unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::Aggregate(double(*fcn)(double)) const
   {
      double res = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr != nullptr) {
            const V& v = itr->GetValue();
            res += fcn(v);
            itr = itr->GetNext();
         }
      }
      return res;
   }

   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   //__host__ __device__ bool AreEqual(unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& a, unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& b)
   //{
   //   if (&a == &b) return true;
   //   return a == b;
   //}

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::iterator(unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& m) : refMap(m), ptr(nullptr), bucket(-1)
   {
      Reset();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::Reset()
   {
      if (refMap.empty()) {
         ptr = nullptr;
         bucket = -1;
         return;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (refMap.buckets[k] != nullptr) {
            bucket = k;
            ptr = refMap.buckets[bucket];
            break;
         }
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::Next()
   {
      if (bucket == -1) {
         return;
      }
      if (ptr != nullptr) {
         ptr = ptr->GetNext();
      }
      if (ptr == nullptr) {
         for (int k = bucket + 1; k < NUM_BUCKETS; ++k) {
            if (refMap.buckets[k] != nullptr) {
               bucket = k;
               ptr = refMap.buckets[bucket];
               return;
            }
         }
         bucket = -1;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::End()
   {
      ptr = nullptr;
      bucket = -1;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::HasNext()
   {
      return bucket != -1;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator++()
   {
      Next();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator==(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator& other) const
   {
      return ptr == other.ptr && bucket == other.bucket;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator!=(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator& other) const
   {
      return !(*this == other);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ LinkedListKV::Node<K, V>& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator*()
   {
      return *ptr;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator=(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator& other)
   {
      refMap = other.refMap;
      ptr = other.ptr;
      bucket = other.bucket;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ K& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::GetKey() const
   {
      if (bucket == -1 || ptr == nullptr) {
         printf("Illegal call to map iterator GetKey(): nullptr pointer, no more element\n");
      }
      return ptr->GetKey();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ V& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::GetValue() const
   {
      if (bucket == -1 || ptr == nullptr) {
         printf("Illegal call to map iterator GetValue(): nullptr pointer, no more element\n");
      }
      return ptr->GetValue();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::begin()
   {
      return unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator(*this);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::end()
   {
      unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator res(*this);
      res.End();
      return res;
   }

   // TODO: implemented properly
   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::cbegin()
   {
      return begin();
   }

   // TODO: implemented properly
   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::cend()
   {
      return end();
   }
}

#endif
