#ifndef STRING_H
#define STRING_H

#include "cuda_runtime.h"

class String
{
public:
   static const int MAX_LEN = 15;

   __host__ __device__ String(char const * s);
   __host__ __device__ ~String();
   __host__ __device__ char* Get();

private:
   char* str;
};

#endif