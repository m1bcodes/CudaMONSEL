//#ifndef _SET_CUH_
//#define _SET_CUH_
//
//#include <cuda_runtime.h>
//
//#include "Hasher.cuh"
//#include "LinkedList.cuh"
//
//namespace Set
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __constant__ const int NUM_BUCKETS = 23;
//#else
//   const int NUM_BUCKETS = 23;
//#endif
//
//   template<typename T>
//   class Set
//   {
//   public:
//      __host__ __device__ Set(Hasher::pHasher);
//      __host__ __device__ void Put(T);
//      __host__ __device__ bool Exists(T, bool cmp(T, T));
//      __host__ __device__ unsigned int Hash(T);
//      __host__ __device__ LinkedList::Node<T>* GetBucket(int n);
//
//   private:
//      LinkedList::Node<T>* buckets[NUM_BUCKETS] = { NULL };
//      Hasher::pHasher hasher;
//   };
//
//   template<typename T>
//   __host__ __device__ Set<T>::Set(Hasher::pHasher hasher) : hasher(hasher)
//   {
//   }
//
//   template<typename T>
//   __host__ __device__ void Set<T>::Put(T v)
//   {
//      auto bucketNum = Hash(v) % NUM_BUCKETS;
//      LinkedList::InsertHead<T>(&buckets[bucketNum], v);
//   }
//
//   template<typename T>
//   __host__ __device__ bool Set<T>::Exists(T v, bool (*cmp)(T, T))
//   {
//      auto bucketIdx = Hash(v) % NUM_BUCKETS;
//      return LinkedList::Exists<T>(buckets[bucketIdx], v, cmp);
//   }
//
//   template<typename T>
//   __host__ __device__ unsigned int Set<T>::Hash(T v)
//   {
//      return hasher((char*)&v, sizeof(v));
//   }
//
//   template<typename T>
//   __host__ __device__ LinkedList::Node<T>* Set<T>::GetBucket(int n)
//   {
//      return buckets[n % NUM_BUCKETS];
//   }
//}
//
//#endif

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
   class Set
   {
   public:
      typedef bool (*pCmp)(T, T);
      __host__ __device__ Set(Hasher::pHasher, pCmp);
      __host__ __device__ void Put(T);
      __host__ __device__ bool Exists(T);
      __host__ __device__ unsigned int Hash(T);
      __host__ __device__ LinkedList::Node<T>* GetBucket(int);
      __host__ __device__ void Remove(T);
      __host__ __device__ void RemoveAll();

   private:
      LinkedList::Node<T>* buckets[NUM_BUCKETS] = { NULL };
      Hasher::pHasher hasher;
      pCmp cmp;
   };

   template<typename T>
   __host__ __device__ Set<T>::Set(Hasher::pHasher hasher, pCmp cmp) : hasher(hasher), cmp(cmp)
   {
   }

   template<typename T>
   __host__ __device__ void Set<T>::Put(T v)
   {
      auto bucketNum = Hash(v) % NUM_BUCKETS;
      if (!Exists(v)) {
         LinkedList::InsertHead<T>(&buckets[bucketNum], v);
      }
   }

   template<typename T>
   __host__ __device__ bool Set<T>::Exists(T v)
   {
      auto bucketIdx = Hash(v) % NUM_BUCKETS;
      return LinkedList::Exists<T>(buckets[bucketIdx], v, cmp);
   }

   template<typename T>
   __host__ __device__ unsigned int Set<T>::Hash(T v)
   {
      return hasher((char*)&v, sizeof(v));
   }

   template<typename T>
   __host__ __device__ LinkedList::Node<T>* Set<T>::GetBucket(int n)
   {
      return buckets[n % NUM_BUCKETS];
   }

   template<typename T>
   __host__ __device__ void Set<T>::Remove(T v)
   {
      //const int badBucket = -1;
      //int bucketIdx = badBucket;
      //for (int k = 0; k < Set::NUM_BUCKETS; ++k) {
      //   auto bucketHead = &buckets[k];
      //   auto bucketItr = bucketHead;
      //   while (bucketItr != NULL) {
      //      auto val = bucketItr->GetValue();
      //      if (cmp(val, v)) {
      //         bucketIdx = k;
      //         break;
      //      }
      //      bucketItr = bucketItr->GetNext();
      //   }
      //   if (bucketIdx != badBucket) {
      //      break;
      //   }
      //}
   }

   template<typename T>
   __host__ __device__ void Set<T>::RemoveAll()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::RemoveAll(&buckets[k]);
      }
   }
}

#endif
