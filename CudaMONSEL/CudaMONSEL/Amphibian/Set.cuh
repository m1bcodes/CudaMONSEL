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

   template<typename T>
   class Iterator;
   
   template<typename T>
   class Set
   {
      friend class Iterator<T>;
   public:
      typedef bool (*pCmp)(T&, T&);
      __host__ __device__ Set(Hasher::pHasher, pCmp);
      __host__ __device__ Set(const Set&);
      __host__ __device__ Set<T>& Set<T>::operator=(const Set<T>&);
      __host__ __device__ ~Set();

      __host__ __device__ void Initialize();
      __host__ __device__ void DeepCopy(const Set&);
      __host__ __device__ void ClearAndCopy(const Set&);
      __host__ __device__ void Put(T);
      __host__ __device__ void Add(Set&);
      __host__ __device__ bool Exists(T);
      __host__ __device__ bool IsEmpty();
      __host__ __device__ int Size();
      __host__ __device__ unsigned int HashCode();
      __host__ __device__ unsigned int Hash(T);
      __host__ __device__ void Remove(T);
      __host__ __device__ void RemoveAll();
      //__host__ __device__ LinkedList::Node<T>* AsList();
      __host__ __device__ LinkedList::Node<T>* GetBucket(int); // DEBUGGING PURPOSES

   private:
      unsigned int GetBucketIdx(T v);
      LinkedList::Node<T>* buckets[NUM_BUCKETS];
      Hasher::pHasher hasher;
      pCmp cmp;
   };

   template<typename T>
   __host__ __device__ Set<T>::Set(Hasher::pHasher hasher, pCmp cmp) : hasher(hasher), cmp(cmp)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
   }

   template<typename T>
   __host__ __device__ Set<T>::Set(const Set& other)
   {
      ClearAndCopy(other);
   }

   template<typename T>
   __host__ __device__ Set<T>& Set<T>::operator=(const Set<T>& other)
   {
      if (this == &other) return;
      ClearAndCopy(other);
      return *this;
   }

   template<typename T>
   __host__ __device__ Set<T>::~Set()
   {
      RemoveAll();
   }

   template<typename T>
   __host__ __device__ void Set<T>::Initialize()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = NULL;
      }
      hasher = NULL;
      cmp = NULL;
   }

   template<typename T>
   __host__ __device__ void Set<T>::DeepCopy(const Set<T>& other)
   {
      RemoveAll();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         auto itr = other.buckets[k];
         while (itr != NULL) {
            Put(itr->GetValue());
            itr = itr->GetNext();
         }
      }
   }

   template<typename T>
   __host__ __device__ void Set<T>::ClearAndCopy(const Set<T>& other)
   {
      Initialize();
      hasher = other.hasher;
      cmp = other.cmp;

      DeepCopy(other);
   }

   template<typename T>
   __host__ __device__ void Set<T>::Put(T v)
   {
      if (!Exists(v)) {
         LinkedList::InsertHead<T>(&buckets[GetBucketIdx(v)], v);
      }
   }

   template<typename T>
   __host__ __device__ void Set<T>::Add(Set<T>& other)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         auto itr = other.buckets[k];
         while (itr != NULL) {
            Put(itr->GetValue());
            itr = itr->GetNext();
         }
      }
   }

   template<typename T>
   __host__ __device__ bool Set<T>::Exists(T v)
   {
      return LinkedList::Exists<T>(buckets[GetBucketIdx(v)], v, cmp);
   }

   template<typename T>
   __host__ __device__ bool Set<T>::IsEmpty()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != NULL) {
            return false;
         }
      }
      return true;
   }

   template<typename T>
   __host__ __device__ int Set<T>::Size()
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

   template<typename T>
   __host__ __device__ unsigned int Set<T>::Hash(T v)
   {
      return hasher((char*)&v, sizeof(v));
   }

   template<typename T>
   __host__ __device__ unsigned int Set<T>::HashCode()
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedList::HashCode(buckets[k], hasher);
      }

      return res;
   }

   template<typename T>
   __host__ __device__ unsigned int Set<T>::GetBucketIdx(T v)
   {
      return Hash(v) % NUM_BUCKETS;
   }

   template<typename T>
   __host__ __device__ LinkedList::Node<T>* Set<T>::GetBucket(int n)
   {
      return buckets[n % NUM_BUCKETS];
   }

   template<typename T>
   __host__ __device__ void Set<T>::Remove(T v)
   {
      LinkedList::Remove(&buckets[GetBucketIdx(v)], v, cmp);
   }

   template<typename T>
   __host__ __device__ void Set<T>::RemoveAll()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::RemoveAll(&buckets[k]);
      }
   }

   //template<typename T>
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

   template<typename T>
   class Iterator
   {
   public:
      __host__ __device__ Iterator(Set<T>&);
      __host__ __device__ void Reset();
      __host__ __device__ void Next();

      __host__ __device__ void operator=(const Iterator&);

      __host__ __device__ bool HasNext();

      __host__ __device__ T GetValue();

   private:
      LinkedList::Node<T>* ptr;
      Set<T>& refSet;
      int bucket;
   };

   template<typename T>
   __host__ __device__ Iterator<T>::Iterator(Set<T>& m) : refSet(m), ptr(NULL), bucket(-1)
   {
      Reset();
   }

   template<typename T>
   __host__ __device__ void Iterator<T>::Reset()
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

   template<typename T>
   __host__ __device__ void Iterator<T>::Next()
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

   template<typename T>
   __host__ __device__ void Iterator<T>::operator=(const Iterator<T>& other)
   {
      ptr = other.ptr;
      bucket = other.bucket;
   }

   template<typename T>
   __host__ __device__ bool Iterator<T>::HasNext()
   {
      return bucket != -1;
   }


   template<typename T>
   __host__ __device__ T Iterator<T>::GetValue()
   {
      return (bucket == -1) ? NULL : ptr->GetValue();
   }
}

#endif
