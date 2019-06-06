#ifndef _UNORDERED_SET_CUH_
#define _UNORDERED_SET_CUH_

#include <cuda_runtime.h>

#include "Hasher.cuh"
#include "LinkedList.cuh"

namespace amp
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static const int NUM_BUCKETS = 23;
#else
   static const int NUM_BUCKETS = 23;
#endif

   template<typename T, typename Hash, typename Pred>
   class Iterator;

   template<typename T, typename Hash, typename Pred>
   class unordered_set
   {
      friend class Iterator<T, Hash, Pred>;
   public:
      __host__ __device__ unordered_set();
      __host__ __device__ unordered_set(const unordered_set&);
      __host__ __device__ unordered_set& operator=(const unordered_set&);
      __host__ __device__ ~unordered_set();

      __host__ __device__ bool operator==(unordered_set&);

      // capacity
      __host__ __device__ bool empty();
      __host__ __device__ int size();

      // iterators

      // element lookup
      __host__ __device__ bool Exists(T&);

      // modifier
      __host__ __device__ void insert(T&);
      __host__ __device__ void erase(T&);
      __host__ __device__ void clear();

      __host__ __device__ void Add(unordered_set&);
      __host__ __device__ unsigned int HashCode();
      __host__ __device__ unsigned int HashCode(T&);
      //__host__ __device__ LinkedList::Node<T>* AsList();

   private:
      __host__ __device__ LinkedList::Node<T>* GetBucket(int);
      __host__ __device__ unsigned int GetBucketIdx(T& v);
      __host__ __device__ void DeepCopy(const unordered_set&);
      __host__ __device__ void ClearAndCopy(const unordered_set&);

      LinkedList::Node<T>* buckets[NUM_BUCKETS];
      Pred cmp;
      Hash hasher;
   };

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unordered_set<T, Hash, Pred>::unordered_set()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = nullptr;
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unordered_set<T, Hash, Pred>::unordered_set(const unordered_set<T, Hash, Pred>& other)
   {
      ClearAndCopy(other);
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unordered_set<T, Hash, Pred>& unordered_set<T, Hash, Pred>::operator=(const unordered_set<T, Hash, Pred>& other)
   {
      if (this == &other) return *this;
      ClearAndCopy(other);
      return *this;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unordered_set<T, Hash, Pred>::~unordered_set()
   {
      clear();
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ bool unordered_set<T, Hash, Pred>::operator==(unordered_set<T, Hash, Pred>& other)
   {
      // TODO: could be faster if just check the buckets or hash, without calling exists, since hash function is standardized
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = buckets[k];
         while (itr != nullptr) {
            T v0 = itr->GetValue();
            if (!other.Exists(v0)) {
               return false;
            }
            itr = itr->GetNext();
         }
      }
      return true;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::DeepCopy(const unordered_set<T, Hash, Pred>& other)
   {
      clear();
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
         LinkedList::Node<T>* itr = other.buckets[k];
         while (itr != nullptr) {
            T v = itr->GetValue();
            insert(v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::ClearAndCopy(const unordered_set<T, Hash, Pred>& other)
   {
      if (&other == this) return;
      // TODO: would cause mem leak if the assigned set is not empty cannot 
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         buckets[k] = nullptr;
      }
      DeepCopy(other);
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::insert(T& v)
   {
      if (!Exists(v)) {
         LinkedList::InsertHead<T>(&buckets[GetBucketIdx(v)], v);
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::Add(unordered_set<T, Hash, Pred>& other)
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = other.buckets[k];
         while (itr != nullptr) {
            T v = itr->GetValue();
            insert(v);
            itr = itr->GetNext();
         }
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ bool unordered_set<T, Hash, Pred>::Exists(T& v)
   {
      return LinkedList::Exists<T, Pred>(buckets[GetBucketIdx(v)], v);
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ bool unordered_set<T, Hash, Pred>::empty()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (buckets[k] != nullptr) {
            return false;
         }
      }
      return true;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ int unordered_set<T, Hash, Pred>::size()
   {
      int c = 0;
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::Node<T>* itr = buckets[k];
         while (itr != nullptr) {
            ++c;
            itr = itr->GetNext();
         }
      }
      return c;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unsigned int unordered_set<T, Hash, Pred>::HashCode(T& v)
   {
      return hasher(v);
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unsigned int unordered_set<T, Hash, Pred>::HashCode()
   {
      unsigned int res = 0;

      for (int k = 0; k < NUM_BUCKETS; ++k) {
         res += LinkedList::HashCode<T, Hash>(buckets[k]);
      }

      return res;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ unsigned int unordered_set<T, Hash, Pred>::GetBucketIdx(T& v)
   {
      return HashCode(v) % NUM_BUCKETS;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ LinkedList::Node<T>* unordered_set<T, Hash, Pred>::GetBucket(int n)
   {
      return buckets[n % NUM_BUCKETS];
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::erase(T& v)
   {
      LinkedList::Remove<T, Pred>(&buckets[GetBucketIdx(v)], v);
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void unordered_set<T, Hash, Pred>::clear()
   {
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         LinkedList::RemoveAll(&buckets[k]);
      }
   }

   //template<typename T, typename Hash, typename Pred>
   //__host__ __device__ LinkedList::Node<T>* unordered_set<T>::AsList()
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

   template<typename T, typename Hash, typename Pred>
   class Iterator
   {
   public:
      __host__ __device__ Iterator(unordered_set<T, Hash, Pred>&);
      __host__ __device__ void Reset();
      __host__ __device__ void Next();

      __host__ __device__ void operator=(const Iterator&);

      __host__ __device__ bool HasNext();

      __host__ __device__ T& GetValue();

   private:
      LinkedList::Node<T>* ptr;
      unordered_set<T, Hash, Pred>& refSet;
      int bucket;
   };

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ Iterator<T, Hash, Pred>::Iterator(unordered_set<T, Hash, Pred>& m) : refSet(m), ptr(NULL), bucket(-1)
   {
      Reset();
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void Iterator<T, Hash, Pred>::Reset()
   {
      if (refSet.empty()) {
         ptr = nullptr;
         bucket = -1;
         return;
      }
      for (int k = 0; k < NUM_BUCKETS; ++k) {
         if (refSet.buckets[k] != nullptr) {
            bucket = k;
            ptr = refSet.buckets[bucket];
            break;
         }
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void Iterator<T, Hash, Pred>::Next()
   {
      if (bucket == -1) {
         return;
      }
      if (ptr != nullptr) {
         ptr = ptr->GetNext();
      }
      if (ptr == nullptr) {
         for (int k = bucket + 1; k < NUM_BUCKETS; ++k) {
            if (refSet.buckets[k] != nullptr) {
               bucket = k;
               ptr = refSet.buckets[bucket];
               return;
            }
         }
         bucket = -1;
      }
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ void Iterator<T, Hash, Pred>::operator=(const Iterator<T, Hash, Pred>& other)
   {
      ptr = other.ptr;
      bucket = other.bucket;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ bool Iterator<T, Hash, Pred>::HasNext()
   {
      return bucket != -1;
   }

   template<typename T, typename Hash, typename Pred>
   __host__ __device__ T& Iterator<T, Hash, Pred>::GetValue()
   {
      if (bucket == -1 || ptr == nullptr) {
         printf("Illegal call to set iterator GetValue(): nullptr, no more element");
      }
      return ptr->GetValue();
   }
}

#endif