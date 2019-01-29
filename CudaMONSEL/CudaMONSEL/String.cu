#include "String.cuh"

#include <stdlib.h>
#include <malloc.h>

__host__ __device__ String::String(char const * s)
{
   for (int k = 0; *s != NULL; ++s, ++k) {
      if (k == MAX_LEN - 1) {
         str[k] = '\0';
         break;
      }
      str[k] = *s;
   }
}

__host__ __device__ char* String::Get()
{
   return str;
}
