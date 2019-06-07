#include "Amphibian/String.cuh"

#include <stdio.h>
#include <string.h>

#include "Amphibian/Hasher.cuh"

namespace amp
{
   __host__ __device__ string::string()
   {
      copy("");
   }

   __host__ __device__ string::string(const string& s)
   {
      copy(s.str);
   }

   __host__ __device__ string::string(char const * s)
   {
      copy(s);
   }

   __host__ __device__ void string::operator=(const string& s)
   {
      if (&s == this) return;
      copy(s.str);
   }

   __host__ __device__ void string::operator=(char const * s)
   {
      copy(s);
   }

   __host__ __device__ bool string::operator==(const string& a) const
   {
      return AreEqual(str, a.str);
   }

   __host__ __device__ const char* string::c_str() const
   {
      return str;
   }

   __host__ __device__ int string::size() const
   {
      int c = 0;
      while (str[c++] != NULL_CHAR);
      return c;
   }

   __host__ __device__ void string::copy(char const * s)
   {
      memset(str, NULL_CHAR, MAX_LEN);
      int k;
      for (k = 0; s[k] != NULL_CHAR; ++k) {
         if (k == MAX_LEN - 1) {
            printf("Length of string exceeded %d.", MAX_LEN);
            break;
         }
         str[k] = s[k];
      }
      str[k] = NULL_CHAR;
   }

   __host__ __device__ unsigned int string::hashcode() const
   {
      return Hasher::Hash(str, size());
   }

   __host__ __device__ void IToA(char * d, int n, int maxArrayLen)
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

   __host__ __device__ int AToI(char const * str)
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

      while (str[idx] != NULL_CHAR) {
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

   __host__ __device__ bool AreEqual(string& a, string& b)
   {
      if (&a == &b) return true;
      return a == b;
   }

   __host__ __device__ bool AreEqual(char const * const a, char const * const b)
   {
      int k = 0;
      while (true) {
         if (a[k] == NULL_CHAR && b[k] == NULL_CHAR) {
            break;
         }
         if (b[k] == NULL_CHAR || a[k] == NULL_CHAR) {
            return false;
         }
         if (a[k] != b[k]) {
            return false;
         }
         ++k;
      }
      return true;
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