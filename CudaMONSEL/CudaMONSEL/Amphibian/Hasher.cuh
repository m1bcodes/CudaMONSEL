/*
**************************************************************************
*                                                                        *
*          General Purpose Hash Function Algorithms Library              *
*                                                                        *
* Author: Arash Partow - 2002                                            *
* URL: http://www.partow.net                                             *
* URL: http://www.partow.net/programming/hashfunctions/index.html        *
*                                                                        *
* Copyright notice:                                                      *
* Free use of the General Purpose Hash Function Algorithms Library is    *
* permitted under the guidelines and in accordance with the MIT License. *
* http://www.opensource.org/licenses/MIT                                 *
*                                                                        *
**************************************************************************
*/

#ifndef _HASHER_CUH_
#define _HASHER_CUH_

#include <cuda_runtime.h>

namespace Hasher
{
   typedef unsigned int(*pHasher)(const char*, unsigned int len);

   __host__ __device__ unsigned int RSHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int JSHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int PJWHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int ELFHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int BKDRHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int SDBMHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int DJBHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int DEKHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int BPHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int FNVHash(const char* str, unsigned int len);
   __host__ __device__ unsigned int APHash(const char* str, unsigned int len);

   __host__ __device__ unsigned int Hash(const char* str, unsigned int len);

   struct IntHashFcn
   {
      __host__ __device__ inline unsigned int operator() (const int& str) const
      {
         return Hash((char *)&str, sizeof(str));
      }
   };

   struct DoubleHashFcn
   {
      __host__ __device__ inline unsigned int operator() (const double& str) const
      {
         return Hash((char *)&str, sizeof(str));
      }
   };
}

namespace Comparator
{
   struct IntCompareFcn
   {
      __host__ __device__ inline bool operator() (const int& lhs, const int& rhs) const
      {
         return lhs == rhs;
      }
   };

   struct DoubleCompareFcn
   {
      __host__ __device__ inline bool operator() (const double& lhs, const double& rhs) const
      {
         return lhs == rhs;
      }
   };
}

#endif
