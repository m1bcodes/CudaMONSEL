#ifndef STRING_CUH
#define STRING_CUH

#include "cuda_runtime.h"

class String
{
public:
   __host__ __device__ static void IToA(char*, int, int maxLen = 11 /* integer limit */);
   __host__ __device__ static bool IsEqual(String, String);


   __host__ __device__ String();
   __host__ __device__ String(char const *);

   __host__ __device__ void operator=(String);

   __host__ __device__ char* Get();

private:
   static const int MAX_LEN = sizeof(char) * 31;

   __host__ __device__ void Copy(char const *);

   char str[MAX_LEN];
};

#endif