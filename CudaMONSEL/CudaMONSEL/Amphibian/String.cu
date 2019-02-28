#include "String.cuh"

#include <stdio.h>

namespace String
{
   __host__ __device__ String::String()
   {
      Copy("");
   }

   __host__ __device__ String::String(String& s)
   {
      Copy(s.Get());
   }

   __host__ __device__ String::String(char const * s)
   {
      Copy(s);
   }

   __host__ __device__ void String::operator=(String s)
   {
      Copy(s.Get());
   }

   __host__ __device__ bool String::operator==(String a)
   {
      int idx = 0;
      while (true) {
         if (str[idx] == NULL && a.str[idx] == NULL) {
            return true;
         }
         if (str[idx] == NULL || a.str[idx] == NULL) {
            return false;
         }
         if (str[idx] != a.str[idx]) {
            return false;
         }
         ++idx;
      }
      return true;
   }

   __host__ __device__ char* String::Get()
   {
      return str;
   }

   __host__ __device__ int String::Length()
   {
      int c = 0;
      while (str[c] != NULL) {
         ++c;
      }
      return c;
   }

   __host__ __device__ void String::Copy(char const * s)
   {
      int k;
      for (k = 0; *s != NULL; ++s, ++k) {
         if (k == MAX_LEN - 1) {
            printf("Length of string exceeded %d.", MAX_LEN);
            break;
         }
         str[k] = *s;
      }
      str[k] = NULL;
   }

   __host__ __device__ void IToA(char* d, int n, int maxArrayLen)
   {
      if (maxArrayLen < 1) {
         return;
      }
      d[0] = NULL;
      if (maxArrayLen == 1) {
         return;
      }

      int idx = 0;
      if (n < 0) {
         d[idx] = '-';
         idx = 1;
         n *= -1;
      }

      do {
         int m = n % 10;
         n /= 10;
         char a = '0' + m;
         d[idx] = a;
         ++idx;
      } while (n != 0 && idx < maxArrayLen);

      if (idx == maxArrayLen) {
         printf("array too small to contain the entire integer\n");
      }
      
      if (d[0] != '-') {
         for (int k = 0; k < idx / 2; ++k) {
            int lastIdx = (idx - 1) - k;
            char a = d[lastIdx];
            d[lastIdx] = d[k];
            d[k] = a;
         }
      }
      else {
         for (int k = 1; k < idx / 2 + 1; ++k) {
            int lastIdx = idx - k;
            char a = d[lastIdx];
            d[lastIdx] = d[k];
            d[k] = a;
         }
      }
      d[idx] = NULL;
   }

   __host__ __device__ int AToI(char* d)
   {
      int mult = 1;
      int idx = 0;
      if (d[0] == '-') {
         mult *= -1;
         idx = 1;
      }

      int res = 0;
      do {
         char di = d[idx];
         if (di < '0' || di > '9') {
            printf("invalid digit");
            return 0;
         }
         int n = di - '0';
         if (res > (MAX_SIGNED_INTEGER - 1 - n) / 10) {
            printf("array contains a number that is out of the integer range\n");
            break;
         }
         res = res * 10 + n;
         ++idx;
      } while (d[idx] != NULL);

      return res*mult;
   }

   __host__ __device__ bool AreEqual(String a, String b)
   {
      return a == b;
   }

   __host__ __device__ bool StartsWith(char* src, char* target)
   {
      int c = 0;
      while (true) {
         if (target[c] == NULL) {
            return true;
         }
         if (src[c] == NULL) {
            return false;
         }
         if (src[c] != target[c]) {
            return false;
         }
         ++c;
      }
   }
}
