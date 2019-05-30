#ifndef _SET_CUH_
#define _SET_CUH_

#include <cuda_runtime.h>

#include "Hasher.cuh"
#include "LinkedList.cuh"

namespace Set
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int NUM_BUCKETS = 23;
#else
   const int NUM_BUCKETS = 23;
#endif

   template<typename T, typename TCompare, typename THasher>
   class Iterator;

   template<typename T, typename TCompare, typename THasher>
   class Set
   {
      friend class Iterator<T, TCompare, THasher>;
   public:
      __host__ __device__ Set();
      __host__ __device__ Set(const Set&);
      __host__ __device__ Set& operator=(const Set&);
      __host__ __device__ ~Set();

      __host__ __device__ bool operator==(Set&);

      __host__ __device__ void Initialize();
      __host__ __device__ void DeepCopy(const Set&);
      __host__ __device__ void ClearAndCopy(const Set&);
      __host__ __device__ void Put(T&);
      __host__ __device__ void Add(Set&);
      __host__ __device__ bool Exists(T&);
      __host__ __device__ bool IsEmpty();
      __host__ __device__ int Size();
      __host__ __device__ unsigned int HashCode();
      __host__ __device__ unsigned int Hash(T&);
      __host__ __device__ void Remove(T&);
      __host__ __device__ void RemoveAll();
      //__host__ __device__ LinkedList::Node<T>* AsList();

   private:
      __host__ __device__ LinkedList::Node<T>* GetBucket(int);
      __host__ __device__ unsigned int GetBucketIdx(T& v);

      LinkedList::Node<T>* buckets[NUM_BUCKETS];
      TCompare cmp;
      THasher hasher;
   };

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ Set<T, TCompare, THasher>::Set()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ Set<T, TCompare, THasher>::Set(const Set<T, TCompare, THasher>& other)
   {
      ClearAndCopy(other);
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ Set<T, TCompare, THasher>& Set<T, TCompare, THasher>::operator=(const Set<T, TCompare, THasher>& other)
   {
      if (this == &other) return *this;
      ClearAndCopy(other);
      return *this;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ Set<T, TCompare, THasher>::~Set()
   {
      RemoveAll();
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ bool Set<T, TCompare, THasher>::operator==(Set<T, TCompare, THasher>& other)
   {
      // TODO: could be faster if just check the buckets or hash, without calling exists, since hash function is standardized
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = buckets[k];
         while (itr != NULL) {
            T v0 = itr->GetValue();
            if (!other.Exists(v0)) {
               return false;
            }
            itr = itr->GetNext();
         }
      }
      return true;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::Initialize()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::DeepCopy(const Set<T, TCompare, THasher>& other)
   {
      RemoveAll();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         LinkedList::Node<T>* itr = other.buckets[k];
         while (itr != NULL) {
            T v = itr->GetValue();
            Put(v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::ClearAndCopy(const Set<T, TCompare, THasher>& other)
   {
      Initialize();

      DeepCopy(other);
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::Put(T& v)
   {
      if (!Exists(v)) {
         LinkedList::InsertHead<T>(&buckets[GetBucketIdx(v)], v);
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::Add(Set<T, TCompare, THasher>& other)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = other.buckets[k];
         while (itr != NULL) {
            T v = itr->GetValue();
            Put(v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ bool Set<T, TCompare, THasher>::Exists(T& v)
   {
      return LinkedList::Exists<T, TCompare>(buckets[GetBucketIdx(v)], v);
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ bool Set<T, TCompare, THasher>::IsEmpty()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != NULL) {
            return false;
         }
      }
      return true;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ int Set<T, TCompare, THasher>::Size()
   {
      int c = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = buckets[k];
         while (itr != NULL) {
            ++c;
            itr = itr->GetNext();
         }
      }
      return c;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ unsigned int Set<T, TCompare, THasher>::Hash(T& v)
   {
      return hasher(v);
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ unsigned int Set<T, TCompare, THasher>::HashCode()
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedList::HashCode<T, THasher>(buckets[k]);
      }

      return res;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ unsigned int Set<T, TCompare, THasher>::GetBucketIdx(T& v)
   {
      return Hash(v) % NUM_BUCKETS;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ LinkedList::Node<T>* Set<T, TCompare, THasher>::GetBucket(int n)
   {
      return buckets[n % NUM_BUCKETS];
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::Remove(T& v)
   {
      LinkedList::Remove<T, TCompare>(&buckets[GetBucketIdx(v)], v);
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Set<T, TCompare, THasher>::RemoveAll()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::RemoveAll(&buckets[k]);
      }
   }

   //template<typename T, typename TCompare, typename THasher>
   //__host__ __device__ LinkedList::Node<T>* Set<T>::AsList()
   //{
   //   LinkedList::Node<T>* res = NULL;
   //   for (int k = 0; k < NUM_BUCKETS; ++k) {
   //      auto itr = buckets[k];
   //      while (itr != NULL) {
   //         LinkedList::InsertHead<T>(&res, itr->GetValue());
   //         itr = itr->GetNext();
   //      }
   //   }
   //   return res;
   //}

   template<typename T, typename TCompare, typename THasher>
   class Iterator
   {
   public:
      __host__ __device__ Iterator(Set<T, TCompare, THasher>&);
      __host__ __device__ void Reset();
      __host__ __device__ void Next();

      __host__ __device__ void operator=(const Iterator&);

      __host__ __device__ bool HasNext();

      __host__ __device__ T& GetValue();

   private:
      LinkedList::Node<T>* ptr;
      Set<T, TCompare, THasher>& refSet;
      int bucket;
   };

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ Iterator<T, TCompare, THasher>::Iterator(Set<T, TCompare, THasher>& m) : refSet(m), ptr(NULL), bucket(-1)
   {
      Reset();
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Iterator<T, TCompare, THasher>::Reset()
   {
      if (refSet.IsEmpty()) {
         ptr = NULL;
         bucket = -1;
         return;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (refSet.buckets[k] != NULL) {
            bucket = k;
            ptr = refSet.buckets[bucket];
            break;
         }
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Iterator<T, TCompare, THasher>::Next()
   {
      if (bucket == -1) {
         return;
      }
      if (ptr != NULL) {
         ptr = ptr->GetNext();
      }
      if (ptr == NULL) {
         for (int k = bucket + 1; k < NUM_BUCKETS; ++k) {
            if (refSet.buckets[k] != NULL) {
               bucket = k;
               ptr = refSet.buckets[bucket];
               return;
            }
         }
         bucket = -1;
      }
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ void Iterator<T, TCompare, THasher>::operator=(const Iterator<T, TCompare, THasher>& other)
   {
      ptr = other.ptr;
      bucket = other.bucket;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ bool Iterator<T, TCompare, THasher>::HasNext()
   {
      return bucket != -1;
   }

   template<typename T, typename TCompare, typename THasher>
   __host__ __device__ T& Iterator<T, TCompare, THasher>::GetValue()
   {
      if (bucket == -1 || ptr == NULL) {
         printf("Illegal call to set iterator GetValue(): NULL pointer, no more element");
      }
      return ptr->GetValue();
   }
}

#endif