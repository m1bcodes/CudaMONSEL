#ifndef _STRING_CUH_
#define _STRING_CUH_

#include "cuda_runtime.h"

class String
{
public:
   __host__ __device__ static void IToA(char*, int, int maxLen = 11 /* integer limit */);

   __host__ __device__ String();
   __host__ __device__ String(char const *);

   __host__ __device__ void operator=(String);

   __host__ __device__ char* Get();

   typedef bool(*pAreEqual)(String, String);
   __host__ __device__ static bool AreEqual(String, String);

private:
   static const int MAX_LEN = sizeof(char) * 31;

   __host__ __device__ void Copy(char const *);

   char str[MAX_LEN];
};

//__global__ void kernel(int n)
//{
//   char n_a[16] = "\0";
//   String::IToA(n_a, n);
//   printf("%s\n", n_a);
//}

#endif