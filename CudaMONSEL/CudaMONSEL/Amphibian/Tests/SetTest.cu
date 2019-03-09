#include "SetTest.cuh"

#include "..\Set.cuh"
#include "..\Hasher.cuh"

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

   __device__ void TestA()
   {
      Hasher::pHasher hasher = Hasher::APHash;

      Set::Set<int> set(hasher, [](int a, int b)
      {
         return a == b;
      });

      for (int k = 0; k < 100; ++k) {
         set.Put(k);
      }
      for (int k = 0; k < 100; ++k) {
         if (!set.Exists(k)) {
            printf("number not found: k\n", k);
         }
      }

      //Print(set);
      printf("SetTest::TestA() completed.\n");

      set.RemoveAll();
   }
}
