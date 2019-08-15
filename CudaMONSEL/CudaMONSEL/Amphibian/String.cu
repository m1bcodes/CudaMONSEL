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
      if (&s == this) return;
      copy(s.str);
   }

   __host__ __device__ string::string(const char s[])
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
      return equal(str, a.str);
   }

   __host__ __device__ bool string::operator!=(const string& a) const
   {
      return !equal(str, a.str);
   }

   __host__ __device__ bool string::operator<(const string& a) const
   {
      int i = 0;
      while (true) {
         if (str[i] == NULL_CHAR && a.str[i] != NULL_CHAR) return true;
         if (str[i] != NULL_CHAR && a.str[i] == NULL_CHAR) return false;
         if (str[i] == NULL_CHAR && a.str[i] == NULL_CHAR) return false;
         if (str[i] < a.str[i]) return true;
         if (str[i] > a.str[i]) return false;
         ++i;
      }
   }

   __host__ __device__ string& string::operator+=(const string& other)
   {
      unsigned int len = size(), otherlen = other.size();
      if (len + otherlen + 1 > MAX_LEN) printf("new string too long\n");
      memcpy(str + len, other.str, sizeof(char) * otherlen);

      return *this;
   }
   
   __host__ __device__ string& string::operator+=(const char* other)
   {
      unsigned int otherlen = 0, len = size();
      while (other[otherlen]) ++otherlen;
      if (otherlen >= MAX_LEN) printf("string too long\n");

      if (len + otherlen + 1 > MAX_LEN) printf("new string too long\n");
      memcpy(str + len, other, sizeof(char) * otherlen);

      return *this;
   }

   __host__ __device__ string& string::operator+=(const char ch)
   {
      unsigned int len = size();
      if (len + 2 > MAX_LEN) printf("new string too long\n");
      memcpy(str + len, &ch, sizeof(char));

      return *this;
   }

   __host__ __device__ const char* string::c_str() const
   {
      return str;
   }

   __host__ __device__ int string::find(const char * target) const
   {
      int len = 0;
      while (target[len]) ++len;

      for (int i = 0; i < MAX_LEN - len; ++i) {
         for (int j = 0; j < len; ++j) {
            if (str[i + j] != target[j]) break;
            if (j == len - 1) return i;
         }
      }

      return npos;
   }

   __host__ __device__ bool string::starts_with(const char * target) const
   {
      unsigned int len = 0;
      while (target[len]) ++len;

      for (int i = 0; i < len; ++i) {
         if (str[i] != target[i] || i == MAX_LEN) return false;
      }
      return true;
   }

   __host__ __device__ string string::substr(size_t pos, size_t len) const
   {
      string newstr;
      len = len == npos ? size() - pos : len;
      memcpy(newstr.str, str + pos, sizeof(char) * len);
      return newstr;
   }

   __host__ __device__ int string::size() const
   {
      int c = 0;
      while (str[c]) ++c;
      return c;
   }

   __host__ __device__ int string::length() const
   {
      return size();
   }

   __host__ __device__ const char& string::at(unsigned int i) const
   {
      return str[i];
   }

   __host__ __device__ char& string::operator[](unsigned int i)
   {
      return str[i];
   }

   __host__ __device__ char* string::data()
   {
      return str;
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

   __host__ __device__ unsigned int string::hashCode() const
   {
      return Hasher::Hash(str, size());
   }

   __host__ __device__ void IToA(int n, char * d, int maxArrayLen)
   {
      if (maxArrayLen < 1) {
         return;
      }
      d[0] = NULL_CHAR;
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

   __host__ __device__ string to_string(int d)
   {
      string res;
      IToA(d, res.data());
      return res;
   }

   __host__ __device__ string to_string(double d)
   {
      return "(double)";
   }

   __host__ __device__ bool equal(const string& a, const string& b)
   {
      if (&a == &b) return true;
      return a == b;
   }

   __host__ __device__ bool equal(char const * const a, char const * const b)
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
         if (target[c] == NULL_CHAR) {
            return true;
         }
         if (src[c] == NULL_CHAR) {
            return false;
         }
         if (src[c] != target[c]) {
            return false;
         }
         ++c;
      }
   }

   __host__ __device__ string operator+(const string& lhs, const string& rhs)
   {
      string newstr(lhs);
      newstr += rhs;
      return newstr;
   }

   __host__ __device__ string operator+(const string& lhs, const char * rhs)
   {
      string newstr(lhs);
      newstr += rhs;
      return newstr;
   }

   __host__ __device__ string operator+(const char * lhs, const string& rhs)
   {
      string newstr(lhs);
      newstr += rhs;
      return newstr;
   }
}