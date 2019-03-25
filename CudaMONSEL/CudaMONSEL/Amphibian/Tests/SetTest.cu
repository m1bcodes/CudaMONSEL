#include "SetTest.cuh"

#include "..\Set.cuh"
#include "..\String.cuh"

#include <stdio.h>

namespace SetTest
{
   __device__ void Print(Set::Set<int> set)
   {
      for (int k = 0; k < Set::NUM_BUCKETS; ++k) {
         auto bucket = set.GetBucket(k);
         while (bucket != NULL) {
            printf("%d, ", bucket->GetValue());
            bucket = bucket->GetNext();
         }
         printf("\n");
      }
   }

   __device__ SetTest::SetTest() : DefaultHasher(Hasher::APHash)
   {
   }

   __device__ void SetTest::TestInt()
   {
      Set::Set<int> set(DefaultHasher, [](int& a, int& b)
      {
         return a == b;
      });

      //for (int k = 0; k < 23; ++k) {
      //   auto b = set.GetBucket(k);
      //   printf("%d b");
      //}

      int maxNum = 100;

      for (int k = 0; k < maxNum; ++k) {
         set.Put(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set.Exists(k)) {
            printf("number not found: k\n", k);
         }
      }

      int num1 = 50, num2 = 100;

      set.Remove(num1);
      set.Remove(num2);

      if (set.Exists(num1)) {
         printf("not removing elements properly: %d\n", num1);
      }
      if (set.Exists(num2)) {
         printf("not removing elements properly: %d\n", num2);
      }

      for (int k = 0; k < maxNum - 20; ++k) {
         set.Remove(k);
         if (set.Exists(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      set.RemoveAll();

      for (int k = 0; k < maxNum; ++k) {
         if (set.Exists(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      printf("SetTest::TestInt() completed.\n");
   }

   __device__ void SetTest::TestInt2()
   {
      Set::Set<int> set1(DefaultHasher, [](int& a, int& b)
      {
         return a == b;
      });

      int maxNum = 100;

      for (int k = 0; k < maxNum; ++k) {
         set1.Put(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set1.Exists(k)) {
            printf("number not found: k\n", k);
         }
      }

      Set::Set<int> set2(DefaultHasher, [](int& a, int& b)
      {
         return a == b;
      });
      for (int k = 0; k < maxNum; ++k) {
         set2.Put(k);
      }
      
      auto h1 = set1.HashCode();
      auto h2 = set2.HashCode();
      if (h1 != h2) {
         printf("HashCodes are different: %u, %u", h1, h2);
      }

      printf("SetTest::TestInt2() completed.\n");
   }

   __device__ void SetTest::TestString()
   {
      Set::Set<String::String> set(DefaultHasher, String::AreEqual);

      String::String a("a");
      String::String b("b");
      String::String a2("a");
      String::String a3("a");
      String::String abc("abc");

      set.Put(a);
      set.Put(b);
      set.Put(abc);

      if (!set.Exists(a2)) {
         printf("does not exist: %s\n", a2.Get());
      }

      if (!set.Exists(a3)) {
         printf("does not exist: %s\n", a2.Get());
      }

      String::String d("d");

      if (set.Exists(d)) {
         printf("exists: %s\n", d.Get());
      }

      Set::Iterator<String::String> itr(set);
      while (itr.HasNext()) {
         //printf("%s ", itr.GetValue().Get());
         itr.Next();
      }
      //printf("\n");

      auto set2 = set;
      set2.Put(String::String("XYZ"));

      if (!(set == set2)) {
         printf("sets are different: %d, %d\n", set.HashCode(), set2.HashCode());
      }

      int c2 = 0;
      Set::Iterator<String::String> itr2(set2);
      while (itr2.HasNext()) {
         //printf("%s ", itr2.GetValue().Get());
         ++c2;
         itr2.Next();
      }
      //printf("\n");
      if (c2 != set2.Size()) {
         printf("different sizes: (%d, %d)", c2, set2.Size());
      }

      printf("SetTest::TestString() completed.\n");
   }
}
