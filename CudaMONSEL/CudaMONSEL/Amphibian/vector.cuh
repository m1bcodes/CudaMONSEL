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

   /*
   * This vector class does not increase size automatically (if full while push_back-ing)
   */
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
      __host__ __device__ vector(int cap, const T&);
      __host__ __device__ vector(const vector&);
      __host__ __device__ void operator=(const vector&);
      __host__ __device__ vector(const T [], const T []);
      __host__ __device__ ~vector();

      // element access
      __host__ __device__ const T& at(const int) const;
      __host__ __device__ T& operator[] (const int);
      __host__ __device__ const T& operator[] (const int) const;
      __host__ __device__ const T* data() const;
      __host__ __device__ T* data();

      //iterator
      __host__ __device__ const_iterator begin() const;
      __host__ __device__ const_iterator end() const;

      // capacity
      __host__ __device__ unsigned int size() const;
      __host__ __device__ bool empty() const;
      __host__ __device__ unsigned int capacity() const;

      // modifiers
      __host__ __device__ void clear(const unsigned int c = VECTOR_INITIAL_SIZE);
      //__host__ __device__ void insert(const int&, const T&);
      __host__ __device__ void erase(const const_iterator&);
      __host__ __device__ void push_back(const T&);
      __host__ __device__ void resize(const unsigned int);
      __host__ __device__ void resize(const unsigned int, const T&);
      __host__ __device__ void assign(const T[], const T[]);
      __host__ __device__ void reserve(const unsigned int);
      __host__ __device__ void set_data(T[], const unsigned int);

      // comparison
      __host__ __device__ bool operator==(const vector<T>&) const;

   private:
      T* vec;
      unsigned int cap;
      unsigned int sz;
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

   // TODO: need to check the reference as well
   template<typename T>
   __host__ __device__ bool vector<T>::const_iterator::operator!=(const const_iterator& other) const
   {
      return i != other.i;
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::const_iterator::operator*() const
   {
      if (i < 0 || i > refvec.size())
         printf("vector<T>::vector::const_iterator::operator*: out of range\n");
      return refvec.vec[i];
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
   __host__ __device__ vector<T>::vector(int cap) : cap(cap), sz(0)
   {
      vec = new T[cap];
   }

   template<typename T>
   __host__ __device__ vector<T>::vector(int cap, const T& v) : cap(cap), sz(0)
   {
      vec = new T[cap];
      for (int i = 0; i < cap; ++i)
         push_back(v); // note: assignment op needed
   }

   template<typename T>
   __host__ __device__ vector<T>::vector(const vector& v) : cap(v.cap), sz(0)
   {
      vec = new T[cap];
      for (int i = 0; i < v.size(); ++i)
         push_back(v.vec[i]);
   }

   template<typename T>
   __host__ __device__ void vector<T>::operator=(const vector& other)
   {
      clear(other.cap);
      for (int i = 0; i < other.sz; ++i)
         push_back(other.vec[i]);
   }

   template<typename T>
   __host__ __device__ vector<T>::vector(const T *start, const T *finish) : cap(max(finish - start, VECTOR_INITIAL_SIZE)), sz(0)
   {
      vec = new T[cap];
      const unsigned int numelem = finish - start;
      for (int i = 0; i < numelem; ++i)
         push_back(*(start + i));
   }

   template<typename T>
   __host__ __device__ vector<T>::~vector()
   {
      delete[] vec;
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::at(const int i) const
   {
      if (!(i >= 0 && i < size())) printf("vector<T>::at: out of range\n");
      return vec[i];
   }

   template<typename T>
   __host__ __device__ T& vector<T>::operator[] (const int i)
   {
      if (!(i >= 0 && i < size())) {
         printf("T& vector<T>::operator[%d]: out of range (%d)\n", i, sz);
         return i < 0 ? vec[0] : vec[sz - 1];
      }
      return vec[i];
   }

   template<typename T>
   __host__ __device__ const T& vector<T>::operator[] (const int i) const
   {
      if (!(i >= 0 && i < size())) {
         printf("const T& vector<T>::operator[%d]: out of range (%d)\n", i, sz);
         return i < 0 ? vec[0] : vec[sz - 1];
      }
      return vec[i];
   }

   template<typename T>
   __host__ __device__ T* vector<T>::data()
   {
      return vec;
   }

   template<typename T>
   __host__ __device__ const T* vector<T>::data() const
   {
      return vec;
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

   template<typename T>
   __host__ __device__ unsigned int vector<T>::size() const
   {
      return sz;
   }

   template<typename T>
   __host__ __device__ bool vector<T>::empty() const
   {
      return sz == 0;
   }

   template<typename T>
   __host__ __device__ unsigned int vector<T>::capacity() const
   {
      return cap;
   }

   template<typename T>
   __host__ __device__ void vector<T>::clear(const unsigned int c)
   {
      delete[] vec;
      cap = c;
      vec = new T[cap];
      sz = 0;
   }

   //template<typename T>
   //__host__ __device__ void vector<T>::insert(const vector<T>::const_iterator& itr, const T& t)
   //{
   //   if (itr >= cap) printf("void vector<T>::insert: out of range");
   //}

   template<typename T>
   __host__ __device__ void vector<T>::erase(const vector<T>::const_iterator& itr)
   {
      if (itr.index() >= sz)
         printf("void vector<T>::erase: out of range\n");
      else {
         const unsigned int idx = itr.index();
         vec[idx].~T();
         for (int i = idx + 1; i < sz; ++i)
            vec[i - 1] = vec[i];
         --sz;
      }
   }

   template<typename T>
   __host__ __device__ void vector<T>::push_back(const T& t)
   {
      if (sz >= cap)
         printf("void vector<T>::push_back: out of range\n");
      else
         vec[sz++] = t; // note: relies on assignment op
   }

   // TODO: if newsz < sz?
   template<typename T>
   __host__ __device__ void vector<T>::resize(const unsigned int newsz)
   {
      //if (cap < newsz) {
      //   T* tmp = new T[newsz];
      //   for (int i = 0; i < sz; ++i) {
      //      tmp[i] = vec[i];
      //   }
      //   clear(newsz);
      //}
      clear(newsz);
      sz = newsz;
   }

   template<typename T>
   __host__ __device__ void vector<T>::resize(const unsigned int newsz, const T& v)
   {
      clear(newsz);
      for (int i = 0; i < newsz; ++i) {
         push_back(v); // note: need assignment op
      }
   }

   template<typename T>
   __host__ __device__ void vector<T>::assign(const T start[], const T end[])
   {
      const unsigned int n = end - start;
      clear(n);

      for (int i = 0; i < n; ++i) {
         push_back(*(start + i));
      }
   }

   template<typename T>
   __host__ __device__ void vector<T>::reserve(const unsigned int newcap)
   {
      clear(newcap);
   }

   template<typename T>
   __host__ __device__ void vector<T>::set_data(T v[], const unsigned int sz)
   {
      delete[] vec;
      this->vec = v;
      this->cap = sz;
      this->sz = sz;
   }

   template<typename T>
   __host__ __device__ bool vector<T>::operator==(const vector<T>& other) const
   {
      if (sz != other.sz) return false;
      for (int i = 0; i < sz; ++i) {
         if (!(vec[i] == other.vec[i])) return false; // note: relies on overloading the operator== for the type T
      }
      return true;
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