#include "String.cu"

#include <stdlib.h>
#include <malloc.h>

String::String(char const * s)
{
   str = (char*)malloc(MAX_LEN + 1);
   int k;
   for (k = 0; *s != NULL; ++s, ++k) {
      if (k > MAX_LEN) {
         break;
      }
      str[k] = *s;
   }
   str[k] = '\0';
}

String::~String()
{
   free(str);
}

char* String::Get()
{
   return str;
}
