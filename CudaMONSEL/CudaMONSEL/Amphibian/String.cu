#include "String.cuh"

__host__ __device__ void String::Copy(char const * s)
{
   int k;
   for (k = 0; *s != NULL; ++s, ++k) {
      if (k == MAX_LEN - 1) {
         break;
      }
      str[k] = *s;
   }
   str[k] = '\0';
}

__host__ __device__ String::String()
{
}

__host__ __device__ String::String(char const * s)
{
   Copy(s);
}

__host__ __device__ void String::operator=(String s)
{
   Copy(s.Get());
}

__host__ __device__ char* String::Get()
{
   return str;
}
