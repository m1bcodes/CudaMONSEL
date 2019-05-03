#include "String.cuh"

#include <stdio.h>
#include <string.h>

namespace String
{
   __host__ __device__ bool AreEqual(char const * const a, char const * const b)
   {
      int k = 0;
      while (true) {
         if (a[k] == NULL && b[k] == NULL) {
            return true;
         }
         if (b[k] == NULL || a[k] == NULL) {
            return false;
         }
         if (a[k] != b[k]) {
            return false;
         }
         ++k;
      }
      return true;
   }

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

   __host__ __device__ void String::operator=(String& s)
   {
      if (&s == this) return;
      Copy(s.Get());
   }

   __host__ __device__ void String::operator=(char const * s)
   {
      Copy(s);
   }

   __host__ __device__ bool String::operator==(const String& a) const
   {
      return AreEqual(str, a.str);
   }

   __host__ __device__ char* String::Get()
   {
      return str;
   }

   __host__ __device__ int String::Size()
   {
      int c = 0;
      while (str[c++] != NULL);
      return c;
   }

   __host__ __device__ void String::Copy(char const * s)
   {
      memset(str, NULL, MAX_LEN);
      int k;
      for (k = 0; s[k] != NULL; ++k) {
         if (k == MAX_LEN - 1) {
            printf("Length of string exceeded %d.", MAX_LEN);
            break;
         }
         str[k] = s[k];
      }
      str[k] = NULL;
   }

   __host__ __device__ unsigned int String::HashCode()
   {
      return Hasher::Hash(str, Size());
   }

   __host__ __device__ void IToA(char* d, int n, int maxArrayLen)
   {
      if (maxArrayLen < 1) {
         return;
      }
      d[0] = '\0';
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
      d[idx] = '\0';
   }

   __host__ __device__ int AToI(char* str)
   {
      int res = 0;
      int mult = 1;

      int idx = 0;
      while (str[idx] == ' ') {
         ++idx;
      }

      if (str[idx] == '-') {
         mult = -1;
         ++idx;
      }
      else if (str[idx] == '+') {
         ++idx;
      }
      else if (!(str[idx] >= '0' && str[idx] <= '9')) {
         printf("invalid digit");
         return 0;
      }

      while (str[idx] != NULL) {
         if (!(str[idx] >= '0' && str[idx] <= '9')) {
            return res;
         }
         int n = (str[idx] - '0') * mult;

         if (mult == 1) {
            if ((INT_MAX - n) / 10 < res) {
               printf("array contains a number that is out of the integer range\n");
               return INT_MAX;
            }
         }
         else {
            if ((INT_MIN - n) / 10 > res) {
               printf("array contains a number that is out of the integer range\n");
               return INT_MIN;
            }
         }
         res = res * 10 + n;

         ++idx;
      }
      return res;
   }

   __host__ __device__ bool AreEqual(String& a, String& b)
   {
      if (&a == &b) return true;
      return a == b;
   }

   __host__ __device__ bool StartsWith(char* src, char* target)
   {
      int c = 0;
      while (true) {
         if (target[c] == '\0') {
            return true;
         }
         if (src[c] == '\0') {
            return false;
         }
         if (src[c] != target[c]) {
            return false;
         }
         ++c;
      }
   }
}
