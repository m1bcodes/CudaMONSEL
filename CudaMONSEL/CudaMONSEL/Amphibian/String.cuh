#ifndef _STRING_CUH_
#define _STRING_CUH_

#include "cuda_runtime.h"

namespace String
{
   class String
   {
   public:
      __host__ __device__ String();
      __host__ __device__ String(char const *);

      __host__ __device__ void operator=(String);
      __host__ __device__ bool operator==(String a);

      __host__ __device__ char* Get();

   private:
      static const int MAX_LEN = sizeof(char) * 31;

      __host__ __device__ void Copy(char const *);

      char str[MAX_LEN];
   };

   __host__ __device__ void IToA(char*, int, int maxLen = 11 /* integer limit */);

   typedef bool(*pStrCmp)(String, String);
   __host__ __device__ bool AreEqual(String, String);

   //__host__ __device__ static void Concatenate(String, String);
}

//__global__ void kernel(int n)
//{
//   char n_a[16] = "\0";
//   String::IToA(n_a, n);
//   printf("%s\n", n_a);
//}

#endif