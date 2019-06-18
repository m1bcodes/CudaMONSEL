#ifndef _VECTOR_CUH_
#define _VECTOR_CUH_

#include <cuda_runtime.h>

#include <stdio.h>

namespace amp
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __device__ static const int VECTOR_INITIAL_SIZE = 16;
#else
   static const int VECTOR_INITIAL_SIZE = 16;
#endif

   template<typename T>
   class vector
   {
   public:
      class const_iterator
      {
      public:
         __host__ __device__ const_iterator(const vector&);
         __host__ __device__ const_iterator(const vector&, const int);
         __host__ __device__ const_iterator(const const_iterator&);

         __host__ __device__ void operator++();
         __host__ __device__ bool operator!=(const const_iterator&) const;
         __host__ __device__ const T& operator*() const;
         __host__ __device__ const T& operator=(const const_iterator&);
         __host__ __device__ bool operator==(const const_iterator&) const;

         __host__ __device__ void begin();
         __host__ __device__ void end();
         __host__ __device__ unsigned int index() const;

      private:
         const vector& refvec;
         int i;
      };

      // *structor
      __host__ __device__ vector(int cap = VECTOR_INITIAL_SIZE);
      __host__ __device__ vector(const vector&);
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
      __host__ __device__ void erase(const const_iterator&);
      __host__ __device__ void push_back(const T&);

      // comparison
      __host__ __device__ bool operator==(const vector<T>&) const;

   private:
      T** vec;
      unsigned int cap;
   };

   template<typename T>
   __host__ __device__ vector<T>::const_iterator::const_iterator(const vector& ref) : refvec(ref), i(0)
   {
   }

   template<typename T>
   __host__ __device__ vector<T>::const_iterator::const_iterator(const vector& ref, const int i) : refvec(ref), i(i)
   {
   }

   template<typename T>
   __host__ __device__ vector<T>::const_iterator::const_iterator(const const_iterator& other) : refvec(other.refvec), i(other.i)
   {  
   }

   template<typename T>
   __host__ __device__ void vector<T>::const_iterator::operator++()
   {
      ++i;
   }

   template<typename T>
   __host__ __device__ bool vector<T>::const_iterator::operator!=(const const_iterator& other) const
   {
      return i != other.i;
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::const_iterator::operator*() const
   {
      if (i < 0 || i > refvec.size()) printf("vector<T>::vector::const_iterator::operator*: out of range\n");
      return *refvec.vec[i];
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::const_iterator::operator=(const const_iterator& other)
   {
      refvec = other.refvec;
      i = other.i;

      return *this;
   }

   template<typename T>
   __host__ __device__ bool vector<T>::const_iterator::operator==(const const_iterator& other) const
   {
      return ((refvec == other.refvec) && (i == other.i));
   }

   template<typename T>
   __host__ __device__ void vector<T>::const_iterator::begin()
   {
      i = 0;
   }

   template<typename T>
   __host__ __device__ void vector<T>::const_iterator::end()
   {
      i = refvec.size();
   }

   template<typename T>
   __host__ __device__ unsigned int vector<T>::const_iterator::index() const
   {
      return i;
   }

   template<typename T>
   __host__ __device__ vector<T>::vector(int cap) : cap(cap)
   {
      vec = new T*[cap];
      for (int i = 0; i < cap; ++i) {
         vec[i] = nullptr;
      }
      //memset(vec, nullptr, sizeof(T*) * cap);
   }

   template<typename T>
   __host__ __device__ vector<T>::vector(const vector& v) : cap(v.cap)
   {
      vec = new T*[cap];
      for (int i = 0; i < cap; ++i) {
         vec[i] = nullptr;
         if (v.vec[i]) push_back(*v.vec[i]);
      }
   }

   template<typename T>
   __host__ __device__ vector<T>::~vector()
   {
      clear();
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::at(const int i) const
   {
      if (!(i >= 0 && i < size())) printf("vector<T>::at out of range\n");
      return *vec[i];
   }

   template<typename T>
   __host__ __device__ T& vector<T>::operator[] (const int i)
   {
      if (!(i >= 0 && i < size())) printf("vector<T>::at out of range\n");
      return *vec[i];
   }

   template<typename T>
   __host__ __device__ vector<T>::const_iterator vector<T>::begin() const
   {
      return vector<T>::const_iterator(*this);
   }

   template<typename T>
   __host__ __device__ vector<T>::const_iterator vector<T>::end() const
   {
      vector<T>::const_iterator res(*this);
      res.end();
      return res;
   }

   // TODO: fix slow
   template<typename T>
   __host__ __device__ unsigned int vector<T>::size() const
   {
      unsigned int i = 0;
      while (vec[i] != nullptr) {
         ++i;
      }
      return i;
   }

   template<typename T>
   __host__ __device__ bool vector<T>::empty() const
   {
      return size() == 0;
   }

   template<typename T>
   __host__ __device__ unsigned int vector<T>::capacity() const
   {
      return cap;
   }

   template<typename T>
   __host__ __device__ void vector<T>::clear()
   {
      unsigned int i = 0;
      while (vec[i] != nullptr) {
         delete vec[i];
         vec[i] = nullptr;
      }
   }

   //template<typename T>
   //__host__ __device__ void vector<T>::insert(const vector<T>::const_iterator& itr, const T& t)
   //{
   //   if (itr >= cap) printf("void vector<T>::insert: out of range");
   //}

   template<typename T>
   __host__ __device__ void vector<T>::erase(const vector<T>::const_iterator& itr)
   {
      unsigned int sz = size();
      if (itr.index() >= sz) printf("void vector<T>::erase: out of range\n");
      else {
         unsigned int idx = itr.index();
         delete vec[idx];
         vec[idx] = nullptr;
         for (int i = idx + 1; i < sz; ++i) {
            vec[i - 1] = vec[i];
         }
         vec[sz - 1] = nullptr;
      }
   }

   template<typename T>
   __host__ __device__ void vector<T>::push_back(const T& t)
   {
      if (size() >= cap) printf("void vector<T>::insert: out of range\n");
      else vec[size()] = new T(t); // note: relies on copy constructor
   }

   template<typename T>
   __host__ __device__ bool vector<T>::operator==(const vector<T>& other) const
   {
      int i = 0;
      while (true) {
         if (i >= cap || i >= other.cap) return false;
         if (!vec[i] && !other.vec[i]) return true;
         if (!vec[i] || !other.vec[i]) return false;
         if (!(*vec[i] == *other.vec[i])) return false; // note: relies on overloading the operator== for the type T
         ++i;
      }
   }

   template<typename T>
   __host__ __device__ typename vector<T>::const_iterator find(const typename vector<T>::const_iterator& start, const typename vector<T>::const_iterator& finish, const T& t)
   {
      typename vector<T>::const_iterator itr(start);
      while (itr != finish) {
         if (*itr == t) return itr; // note: relies on overloading the operator== for the type T
         ++itr;
      }
      itr.end();
      return itr;
   }
}

#endif