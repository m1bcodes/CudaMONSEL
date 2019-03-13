/**
 * An immutable string class.
 */
#ifndef _STRING_CUH_
#define _STRING_CUH_

#include <cuda_runtime.h>
#include "LinkedList.cuh"

namespace String
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int MAX_LEN = sizeof(char) * 32;
#else
   const int MAX_LEN = sizeof(char) * 32;
#endif

   class String
   {
   public:
      __host__ __device__ String();
      __host__ __device__ String(String&);
      __host__ __device__ String(char const *);

      __host__ __device__ void operator=(String);
      __host__ __device__ bool operator==(String);

      __host__ __device__ char* Get();
      __host__ __device__ int Length();

   private:
      __host__ __device__ void Copy(char const *);

      char str[MAX_LEN];
   };

   __host__ __device__ void IToA(char*, int, int maxArrayLen = 11 /* integer limit */);
   __host__ __device__ int AToI(char*);
   //__host__ __device__ float AToF(char*);
   //__host__ __device__ double AToD(char*);
   template<typename T>
   __host__ __device__ T AToF(char* d)
   {
      int mult = 1;
      int idx = 0;
      if (d[0] == '-') {
         mult *= -1;
         idx = 1;
      }

      T res = 0;
      T decDivs = 10;
      bool foundDec = false;
      do {
         char di = d[idx];
         if (di == '.') {
            foundDec = true;
         }
         else if (di >= '0' || di <= '9') {
            int n = di - '0';
            if (foundDec) {
               res = res + n / decDivs;
               decDivs *= 10;
            }
            else {
               res = res * 10 + n;
            }
         }
         ++idx;
      } while (d[idx] != NULL);

      return res*mult;
   }

   typedef bool(*pStrCmp)(String, String);
   __host__ __device__ bool AreEqual(String, String);
   __host__ __device__ bool StartsWith(char* src, char* target);

   //__host__ __device__ static void Concatenate(String, String);
}

//__global__ void kernel(int n)
//{
//   char n_a[16] = "\0";
//   String::IToA(n_a, n);
//   printf("%s\n", n_a);
//}

#endif