#ifndef _UNORDERED_MAP_CUH_
#define _UNORDERED_MAP_CUH_

#include "Amphibian/unordered_set.cuh"
#include "Amphibian/LinkedList.cuh"
#include "Amphibian/Hasher.cuh"

#include <cuda_runtime.h>

#include <stdio.h>

namespace amp
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static const int NUM_MAP_BUCKETS = 23;
#else
   static const int NUM_MAP_BUCKETS = 23;
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

      // element access
      __host__ __device__ V& operator[](const K&);
      __host__ __device__ const V& operator[](const K&) const;

      // modifiers
      __host__ __device__ void insert(LinkedListKVPtr);
      __host__ __device__ bool erase(const K&, V&);
      __host__ __device__ bool erase(const K&);
      __host__ __device__ void clear();

      // hash
      __host__ __device__ unsigned int hashCode() const;
      __host__ __device__ KHasher hash_function() const; // key has function, as per STL

      __host__ __device__ void Initialize();
      __host__ __device__ void ClearAndCopy(const unordered_map&);
      __host__ __device__ void Put(const K&, const V&);
      __host__ __device__ bool ContainsKey(const K&) const;
      __host__ __device__ bool GetValue(const K&, V&) const;
      __host__ __device__ amp::unordered_set<K, KCompare, KHasher> GetKeys();
      __host__ __device__ unsigned int Hash(const K&) const;
      __host__ __device__ void DeepCopy(const unordered_map& other);
      __host__ __device__ double Aggregate(double(*fcn)(double)) const;

      class iterator
      {
      public:
         __host__ __device__ iterator(unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>&);
         __host__ __device__ void Reset();
         __host__ __device__ void Next();
         __host__ __device__ void End();

         __host__ __device__ bool HasNext() const;

         __host__ __device__ void operator++();
         __host__ __device__ bool operator!=(const iterator&) const;
         __host__ __device__ bool operator==(const iterator&) const;
         __host__ __device__ void operator=(const iterator&);

         __host__ __device__ K& GetKey() const;
         __host__ __device__ V& GetValue() const;

         __host__ __device__ LinkedListKV::Node<K, V>& operator*();
         __host__ __device__ LinkedListKV::Node<K, V>* operator->();

      protected:
         unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& refMap;
         LinkedListKV::Node<K, V>* ptr;
         int bucket;
      };

      class const_iterator
      {
      public:
         __host__ __device__ const_iterator(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>&);
         __host__ __device__ void Reset();
         __host__ __device__ void Next();
         __host__ __device__ void End();

         __host__ __device__ bool HasNext() const;

         __host__ __device__ void operator++();
         __host__ __device__ bool operator!=(const const_iterator&) const;
         __host__ __device__ bool operator==(const const_iterator&) const;
         __host__ __device__ void operator=(const const_iterator&);

         __host__ __device__ const K& GetKey() const;
         __host__ __device__ const V& GetValue() const;

         __host__ __device__ const LinkedListKV::Node<K, V>& operator*() const;
         __host__ __device__ const LinkedListKV::Node<K, V>* operator->() const;

      protected:
         const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& refMap;
         const LinkedListKV::Node<K, V>* ptr;
         int bucket;
      };

      // element lookup
      __host__ __device__ iterator find(const K& k);
      __host__ __device__ const_iterator find(const K& k) const;

      // iterators
      __host__ __device__ iterator begin();
      __host__ __device__ iterator end();
      __host__ __device__ const_iterator begin() const;
      __host__ __device__ const_iterator end() const;
      __host__ __device__ const_iterator cbegin() const;
      __host__ __device__ const_iterator cend() const;

      //class Wrapper
      //{
      //public:
      //   __host__ __device__ First(unordered_map& refMap) : refMap(refMap) {}
      //   __host__ __device__ operator V&() const { return (KeyT&)node.key; }
      //private:
      //   unordered_map& refMap;
      //} first;

   private:
      __host__ __device__ int unsigned GetBucketIdx(const K& k) const;
      __host__ __device__ LinkedListKVPtr GetBucket(int); // DEBUGGING PURPOSES

      LinkedListKVPtr buckets[NUM_MAP_BUCKETS];
      V defaultv; // relies on default constrtor
   };

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::unordered_map()
   {
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
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
      clear();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::operator==(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& other) const
   {
      VCompare vcmp;
      if (this == &other) return true;

      if (size() != other.size()) {
         return false;
      }
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         LinkedListKVPtr itr = other.buckets[k];
         while (itr) {
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
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         buckets[k] = nullptr;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::DeepCopy(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
   {
      clear();
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         LinkedListKVPtr itr = other.buckets[k];
         while (itr) {
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
      return LinkedListKV::GetValue<K, V, KCompare>(buckets[GetBucketIdx(k)], k, v);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ amp::unordered_set<K, KCompare, KHasher> unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeys()
   {
      amp::unordered_set<K, KCompare, KHasher> res;
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr) {
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
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::clear()
   {
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
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
      auto itr = find(node->GetKey());
      if (itr == end()) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(node->GetKey())], node);
      }
      else {
         (*itr).MapVal(node->GetValue(), GetFirstParam<V>);
         delete node;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::Put(const K& k, const V& v)
   {
      auto itr = find(k);
      if (itr == end()) {
         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(k)], k, v);
      }
      else {
         (*itr).MapVal(v, GetFirstParam<V>);
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::ContainsKey(const K& k) const
   {
      return LinkedListKV::ContainsKey<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::find(const K& k)
   {
      KCompare kcmp;
      iterator itr(*this);
      while (itr != end()) {
         if (kcmp(itr.GetKey(), k)) return itr;
         ++itr;
      }
      return itr;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::find(const K& k) const
   {
      KCompare kcmp;
      const_iterator itr(*this);
      while (itr != end()) {
         if (kcmp(itr.GetKey(), k)) return itr;
         ++itr;
      }
      return itr;
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
      return Hash(v) % NUM_MAP_BUCKETS;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::hashCode() const
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         res += LinkedListKV::HashCode<K, V, KHasher, VHasher>(buckets[k]);
      }

      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::empty() const
   {
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         if (buckets[k]) {
            return false;
         }
      }
      return true;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ int unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::size() const
   {
      int c = 0;
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr) {
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
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         LinkedListKVPtr itr = buckets[k];
         while (itr) {
            const V& v = itr->GetValue();
            res += fcn(v);
            itr = itr->GetNext();
         }
      }
      return res;
   }   
   
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
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         if (refMap.buckets[k]) {
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
      if (ptr) {
         ptr = ptr->GetNext();
      }
      if (ptr == nullptr) {
         for (int k = bucket + 1; k < NUM_MAP_BUCKETS; ++k) {
            if (refMap.buckets[k]) {
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
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::HasNext() const
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
   __host__ __device__ LinkedListKV::Node<K, V>* unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::iterator::operator->()
   {
      return ptr;
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

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::const_iterator(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>& m) : refMap(m), ptr(nullptr), bucket(-1)
   {
      Reset();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::Reset()
   {
      if (refMap.empty()) {
         ptr = nullptr;
         bucket = -1;
         return;
      }
      for (int k = 0; k < NUM_MAP_BUCKETS; ++k) {
         if (refMap.buckets[k]) {
            bucket = k;
            ptr = refMap.buckets[bucket];
            break;
         }
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::Next()
   {
      if (bucket == -1) {
         return;
      }
      if (ptr) {
         ptr = ptr->GetNext();
      }
      if (ptr == nullptr) {
         for (int k = bucket + 1; k < NUM_MAP_BUCKETS; ++k) {
            if (refMap.buckets[k]) {
               bucket = k;
               ptr = refMap.buckets[bucket];
               return;
            }
         }
         bucket = -1;
      }
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::End()
   {
      ptr = nullptr;
      bucket = -1;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::HasNext() const
   {
      return bucket != -1;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator++()
   {
      Next();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator==(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator& other) const
   {
      return ptr == other.ptr && bucket == other.bucket;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ bool unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator!=(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator& other) const
   {
      return !(*this == other);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ void unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator=(const unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator& other)
   {
      refMap = other.refMap;
      ptr = other.ptr;
      bucket = other.bucket;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ const LinkedListKV::Node<K, V>& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator*() const
   {
      return *ptr;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ const LinkedListKV::Node<K, V>* unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::operator->() const
   {
      return ptr;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ const K& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::GetKey() const
   {
      if (bucket == -1 || ptr == nullptr) {
         printf("Illegal call to map iterator GetKey(): nullptr pointer, no more element\n");
      }
      return ptr->GetKey();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ const V& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator::GetValue() const
   {
      if (bucket == -1 || ptr == nullptr) {
         printf("Illegal call to map iterator GetValue(): nullptr pointer, no more element\n");
      }
      return ptr->GetValue();
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::begin() const
   {
      return unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator(*this);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::end() const
   {
      unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator res(*this);
      res.End();
      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::cbegin() const
   {
      return unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator(*this);
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::cend() const
   {
      unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::const_iterator res(*this);
      res.End();
      return res;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ KHasher unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::hash_function() const
   {
      KHasher h;
      return h;
   }

   // TODO: fix map[key] = val not working
   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ V& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::operator[](const K& k)
   {
      KCompare kcmp;
      for (auto p : *this) {
         if (kcmp(p.first, k)) {
            return p.second;
         }
      }
      return defaultv;
   }

   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
   __host__ __device__ const V& unordered_map<K, V, KCompare, VCompare, KHasher, VHasher>::operator[](const K& k) const
   {
      return (*this)[k];
   }
}

#endif
