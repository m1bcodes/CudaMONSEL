//#ifndef _UNSORTED_MAP_CUH_
//#define _UNSORTED_MAP_CUH_
//
//#include "Amphibian/Set.cuh"
//#include "Amphibian/LinkedList.cuh"
//#include "Amphibian/Hasher.cuh"
//
//#include <cuda_runtime.h>
//
//#include <stdio.h>
//
//namespace Map
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __constant__ static const int NUM_BUCKETS = 23;
//#else
//   static const int NUM_BUCKETS = 23;
//#endif
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   class Iterator;
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   class Map
//   {
//      friend class Iterator<K, V, KCompare, VCompare, KHasher, VHasher>;
//      typedef LinkedListKV::Node<K, V>* LinkedListKVPtr;
//
//   public:
//      __host__ __device__ Map();
//      __host__ __device__ Map(const Map&);
//      __host__ __device__ Map(int);
//      __host__ __device__ Map& operator=(const Map&);
//      __host__ __device__ ~Map();
//
//      __host__ __device__ bool operator==(Map&);
//
//      __host__ __device__ void Initialize();
//      __host__ __device__ void ClearAndCopy(const Map&);
//      __host__ __device__ void Put(K&, V&);
//      __host__ __device__ bool ContainsKey(K&);
//      __host__ __device__ bool GetValue(K&, V&);
//      __host__ __device__ Set::Set<K, KCompare, KHasher> GetKeys();
//      __host__ __device__ unsigned int Hash(K&);
//      __host__ __device__ unsigned int HashCode();
//      __host__ __device__ void DeepCopy(const Map& other);
//      __host__ __device__ bool Remove(K&, V&);
//      __host__ __device__ bool Remove(K&);
//      __host__ __device__ void RemoveAll();
//      __host__ __device__ bool IsEmpty();
//      __host__ __device__ int Size();
//      __host__ __device__ double Aggregate(double(*fcn)(double));
//      //__host__ __device__ LinkedListKV::Node<K, V>* AsList();
//
//      //__host__ __device__ Hasher::pHasher GetHasher(); // DEBUGGING PURPOSES
//      //__host__ __device__ pKeyCmp GetKeyCmp(); // DEBUGGING PURPOSES
//      //__host__ __device__ pValCmp GetValCmp(); // DEBUGGING PURPOSES
//
//   private:
//      __host__ __device__ int unsigned GetBucketIdx(K k);
//      __host__ __device__ LinkedListKVPtr GetBucket(int); // DEBUGGING PURPOSES
//
//      LinkedListKVPtr buckets[NUM_BUCKETS];
//      KCompare kcmp;
//      VCompare vcmp;
//      KHasher khasher;
//      VHasher vhasher;
//   };
//
//   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   //__host__ __device__ Hasher::pHasher Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetHasher()
//   //{
//   //   return hasher;
//   //}
//
//   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   //__host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::pKeyCmp Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeyCmp()
//   //{
//   //   return kcmp;
//   //}
//
//   //template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   //__host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::pValCmp Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValCmp()
//   //{
//   //   return vcmp;
//   //}
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map()
//   {
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         buckets[k] = NULL;
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
//   {
//      //printf("called cc\n");
//      ClearAndCopy(m);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::Map(int a)
//   {
//      Initialize();
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>& Map<K, V, KCompare, VCompare, KHasher, VHasher>::operator=(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& m)
//   {
//      //printf("called =\n");
//      if (&m == this) {
//         return *this;
//      }
//      ClearAndCopy(m);
//      return *this;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::~Map()
//   {
//      //printf("called des\n");
//      RemoveAll();
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::operator==(Map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
//   {
//      if (this == &other) return true;
//
//      if (Size() != other.Size()) {
//         return false;
//      }
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         LinkedListKVPtr itr = other.buckets[k];
//         while (itr != NULL) {
//            K k1 = itr->GetKey();
//            V v0, v1;
//            V found0 = GetValue(k1, v0);
//            V found1 = other.GetValue(k1, v1);
//            if (found0 && !found1 || !found0 && found1) { // only one of the keys is NULL
//               return false;
//            }
//            if (found0 && found1) { // both values are not NULL
//               if (!vcmp(v0, v1)) { // values are different
//                  return false;
//               }
//            }
//            itr = itr->GetNext();
//         }
//      }
//      return true;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Initialize()
//   {
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         buckets[k] = NULL;
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::DeepCopy(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
//   {
//      RemoveAll();
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         // LinkedListKV::DeepCopy(&buckets, other.buckets[k]); hash function may be different
//         LinkedListKVPtr itr = other.buckets[k];
//         while (itr != NULL) {
//            K k = itr->GetKey();
//            V v = itr->GetValue();
//            Put(k, v);
//            itr = itr->GetNext();
//         }
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::ClearAndCopy(const Map<K, V, KCompare, VCompare, KHasher, VHasher>& other)
//   {
//      Initialize();
//
//      DeepCopy(other);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetValue(K& k, V& v)
//   {
//      //LinkedListKV::GetValue<K, V>(buckets[0], k, kcmp);
//      return LinkedListKV::GetValue<K, V, KCompare>(buckets[GetBucketIdx(k)], k, v);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Set::Set<K, KCompare, KHasher> Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetKeys()
//   {
//      Set::Set<K, KCompare, KHasher> res;
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         LinkedListKVPtr itr = buckets[k];
//         while (itr != NULL) {
//            res.Put(itr->GetKey());
//            itr = itr->GetNext();
//         }
//      }
//      return res;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::Remove(K& k, V& ret)
//   {
//      return LinkedListKV::Remove<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, ret);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::Remove(K& k)
//   {
//      V v;
//      return Remove(k, v);
//      //return LinkedListKV::Remove<K, V, KCompare>(&buckets[GetBucketIdx(k)], k, v);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::RemoveAll()
//   {
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         LinkedListKV::RemoveAll<K, V>(&buckets[k]);
//      }
//   }
//
//   template<typename V>
//   __host__ __device__ V GetFirstParam(V& a, V& b) {
//      return a;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Map<K, V, KCompare, VCompare, KHasher, VHasher>::Put(K& k, V& v)
//   {
//      if (!ContainsKey(k)) {
//         LinkedListKV::InsertHead<K, V>(&buckets[GetBucketIdx(k)], k, v);
//      }
//      else {
//         LinkedListKVPtr bucketItr = buckets[GetBucketIdx(k)];
//         while (bucketItr != NULL) {
//            K tmpK = bucketItr->GetKey();
//            if (kcmp(tmpK, k)) {
//               bucketItr->MapVal(v, GetFirstParam<V>);
//               break;
//            }
//            bucketItr = bucketItr->GetNext();
//         }
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::ContainsKey(K& k)
//   {
//      return LinkedListKV::ContainsKey<K, V, KCompare>(buckets[GetBucketIdx(k)], k);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::Hash(K& k)
//   {
//      return khasher(k);
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucketIdx(K v)
//   {
//      return Hash(v) % NUM_BUCKETS;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ unsigned int Map<K, V, KCompare, VCompare, KHasher, VHasher>::HashCode()
//   {
//      unsigned int res = 0;
//
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         res += LinkedListKV::HashCode<K, V, KHasher, VHasher>(buckets[k]);
//      }
//
//      return res;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Map<K, V, KCompare, VCompare, KHasher, VHasher>::IsEmpty()
//   {
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         if (buckets[k] != NULL) {
//            return false;
//         }
//      }
//      return true;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ int Map<K, V, KCompare, VCompare, KHasher, VHasher>::Size()
//   {
//      int c = 0;
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         LinkedListKVPtr itr = buckets[k];
//         while (itr != NULL) {
//            ++c;
//            itr = itr->GetNext();
//         }
//      }
//      return c;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Map<K, V, KCompare, VCompare, KHasher, VHasher>::LinkedListKVPtr Map<K, V, KCompare, VCompare, KHasher, VHasher>::GetBucket(int n)
//   {
//      return buckets[n];
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ double Map<K, V, KCompare, VCompare, KHasher, VHasher>::Aggregate(double(*fcn)(double))
//   {
//      double res = 0;
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         LinkedListKVPtr itr = buckets[k];
//         while (itr != NULL) {
//            V v = itr->GetValue();
//            res += fcn(v);
//            itr = itr->GetNext();
//         }
//      }
//      return res;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool AreEqual(Map<K, V, KCompare, VCompare, KHasher, VHasher>& a, Map<K, V, KCompare, VCompare, KHasher, VHasher>& b)
//   {
//      if (&a == &b) return true;
//      return a == b;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   class Iterator
//   {
//   public:
//      __host__ __device__ Iterator(Map<K, V, KCompare, VCompare, KHasher, VHasher>&);
//      __host__ __device__ void Reset();
//      __host__ __device__ void Next();
//
//      __host__ __device__ bool HasNext();
//
//      __host__ __device__ K GetKey();
//      __host__ __device__ V GetValue();
//
//   private:
//      LinkedListKV::Node<K, V>* ptr;
//      Map<K, V, KCompare, VCompare, KHasher, VHasher>& refMap;
//      int bucket;
//   };
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::Iterator(Map<K, V, KCompare, VCompare, KHasher, VHasher>& m) : refMap(m), ptr(NULL), bucket(-1)
//   {
//      Reset();
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::Reset()
//   {
//      if (refMap.IsEmpty()) {
//         ptr = NULL;
//         bucket = -1;
//         return;
//      }
//      for (int k = 0; k < NUM_BUCKETS; ++k) {
//         if (refMap.buckets[k] != NULL) {
//            bucket = k;
//            ptr = refMap.buckets[bucket];
//            break;
//         }
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ void Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::Next()
//   {
//      if (bucket == -1) {
//         return;
//      }
//      if (ptr != NULL) {
//         ptr = ptr->GetNext();
//      }
//      if (ptr == NULL) {
//         for (int k = bucket + 1; k < NUM_BUCKETS; ++k) {
//            if (refMap.buckets[k] != NULL) {
//               bucket = k;
//               ptr = refMap.buckets[bucket];
//               return;
//            }
//         }
//         bucket = -1;
//      }
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ bool Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::HasNext()
//   {
//      return bucket != -1;
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ K Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::GetKey()
//   {
//      if (bucket == -1 || ptr == NULL) {
//         printf("Illegal call to map iterator GetKey(): NULL pointer, no more element");
//      }
//      return ptr->GetKey();
//   }
//
//   template<typename K, typename V, typename KCompare, typename VCompare, typename KHasher, typename VHasher>
//   __host__ __device__ V Iterator<K, V, KCompare, VCompare, KHasher, VHasher>::GetValue()
//   {
//      if (bucket == -1 || ptr == NULL) {
//         printf("Illegal call to map iterator GetValue(): NULL pointer, no more element");
//      }
//      return ptr->GetValue();
//   }
//}
//
//#endif
