#ifndef _VECTOR_CUH_
#define _VECTOR_CUH_

#include <cuda_runtime.h>

#include <stdio.h>

namespace amp
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static const int VECTOR_INITIAL_SIZE = 23;
#else
   static const int VECTOR_INITIAL_SIZE = 16;
#endif

   template<typename T, typename TCompare>
   class vector
   {
   public:
      class const_iterator
      {
      public:
         __host__ __device__ const_iterator(const vector<T, TCompare>&);
         __host__ __device__ const_iterator(const vector<T, TCompare>&, const int);
         __host__ __device__ const_iterator(const const_iterator&);

         __host__ __device__ void operator++();
         __host__ __device__ bool operator!=(const const_iterator&) const;
         __host__ __device__ const T& operator*() const;

         __host__ __device__ void begin();
         __host__ __device__ void end();

      private:
         const vector& refvec;
         int i;
      };

      // *structor
      __host__ __device__ vector(int cap = VECTOR_INITIAL_SIZE);
      __host__ __device__ ~vector();

      // element access
      __host__ __device__ const T& at(const int) const;
      __host__ __device__ T& operator[] (const int);

      //iterator
      __host__ __device__ const_iterator begin() const;
      __host__ __device__ const_iterator end() const;

      // capacity
      __host__ __device__ unsigned int size() const;
      __host__ __device__ bool empty() const;
      __host__ __device__ unsigned int capacity() const;

      // modifiers
      __host__ __device__ void clear();
      //__host__ __device__ void insert(const int&, const T&);
      //__host__ __device__ void erase(const int&);
      __host__ __device__ void push_back(const T&);

   private:
      T** vec;
      unsigned int cap;
   };

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::const_iterator::const_iterator(const vector& ref) : refvec(ref), i(0)
   {
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::const_iterator::const_iterator(const vector& ref, const int i) : refvec(ref), i(i)
   {
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::const_iterator::const_iterator(const const_iterator& other) : refvec(other.refvec), i(other.i)
   {
   }

   template<typename T, typename TCompare>
   __host__ __device__ void vector<T, TCompare>::const_iterator::operator++()
   {
      ++i;
   }

   template<typename T, typename TCompare>
   __host__ __device__ bool vector<T, TCompare>::const_iterator::operator!=(const const_iterator& other) const
   {
      return i != other.i;
   }

   template<typename T, typename TCompare>
   __host__ __device__ const T& vector<T, TCompare>::const_iterator::operator*() const
   {
      if (i < 0 || i > refvec.size()) printf("vector<T, TCompare>::vector::const_iterator::operator*: out of range\n");
      return *refvec.vec[i];
   }

   template<typename T, typename TCompare>
   __host__ __device__ void vector<T, TCompare>::const_iterator::begin()
   {
      i = 0;
   }

   template<typename T, typename TCompare>
   __host__ __device__ void vector<T, TCompare>::const_iterator::end()
   {
      i = refvec.size();
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::vector(int cap) : cap(cap)
   {
      vec = new T*[cap];
      for (int i = 0; i < cap; ++i) {
         vec[i] = nullptr;
      }
      //memset(vec, nullptr, sizeof(T*) * cap);
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::~vector()
   {
      clear();
   }

   template<typename T, typename TCompare>
   __host__ __device__ const T& vector<T, TCompare>::at(const int i) const
   {
      if (!(i >= 0 && i < size())) printf("vector<T, TCompare>::at out of range\n");
      return *vec[i];
   }

   template<typename T, typename TCompare>
   __host__ __device__ T& vector<T, TCompare>::operator[] (const int i)
   {
      if (!(i >= 0 && i < size())) printf("vector<T, TCompare>::at out of range\n");
      return *vec[i];
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::const_iterator vector<T, TCompare>::begin() const
   {
      return vector<T, TCompare>::const_iterator(*this);
   }

   template<typename T, typename TCompare>
   __host__ __device__ vector<T, TCompare>::const_iterator vector<T, TCompare>::end() const
   {
      vector<T, TCompare>::const_iterator res(*this);
      res.end();
      return res;
   }

   // TODO: fix slow
   template<typename T, typename TCompare>
   __host__ __device__ unsigned int vector<T, TCompare>::size() const
   {
      unsigned int i = 0;
      while (vec[i] != nullptr) {
         ++i;
      }
      return i;
   }

   template<typename T, typename TCompare>
   __host__ __device__ bool vector<T, TCompare>::empty() const
   {
      return size() == 0;
   }

   template<typename T, typename TCompare>
   __host__ __device__ unsigned int vector<T, TCompare>::capacity() const
   {
      return cap;
   }

   template<typename T, typename TCompare>
   __host__ __device__ void vector<T, TCompare>::clear()
   {
      unsigned int i = 0;
      while (vec[i] != nullptr) {
         delete vec[i];
         vec[i] = nullptr;
      }
   }

   //template<typename T, typename TCompare>
   //__host__ __device__ void vector<T, TCompare>::insert(const vector<T, TCompare>::const_iterator& itr, const T& t)
   //{
   //   if (itr >= cap) printf("void vector<T, TCompare>::insert: out of range");
   //}

   template<typename T, typename TCompare>
   __host__ __device__ void vector<T, TCompare>::push_back(const T& t)
   {
      if (size() > cap) printf("void vector<T, TCompare>::insert: out of range\n");
      else vec[size()] = new T(t); // note: relies on copy constructor
   }
}

#endif