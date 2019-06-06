#include "Amphibian/Tests/SetTest.cuh"

#include "Amphibian/unordered_set.cuh"
#include "Amphibian/String.cuh"

#include <stdio.h>

namespace SetTest
{
   struct IntCompare
   {
      __device__ inline bool operator() (const int& lhs, const int& rhs)
      {
         if (&rhs == &lhs) return true;
         return lhs == rhs;
      }
   };

   //__device__ void Print(Set::Set<int, IntCompare> set)
   //{
   //   for (int k = 0; k < Set::NUM_BUCKETS; ++k) {
   //      auto bucket = set.GetBucket(k);
   //      while (bucket != NULL) {
   //         printf("%d, ", bucket->GetValue());
   //         bucket = bucket->GetNext();
   //      }
   //      printf("\n");
   //   }
   //}

   __device__ SetTest::SetTest()
   {
   }

   __device__ void SetTest::TestIntBasic()
   {
      int a = 1, b = 1;
      IntCompare cmp;
      if (!cmp(a, b)) {
         printf("bad: %d, %d\n", a, b);
      }
      amp::unordered_set<int, Hasher::IntHashFcn, IntCompare> set;
      int k = 0;
      set.insert(k);
      set.insert(k);
      for (int k = 0; k < 100; ++k) {
         set.insert(k);
      }

      //Set::Iterator<int, IntCompare> itr(set);
      //while (itr.HasNext()) {
      //   printf("%d ", itr.GetValue());
      //   itr.Next();
      //}
      //printf("\n");

      printf("SetTest::TestIntBasic() completed.\n");
   }

   __device__ void SetTest::TestInt()
   {
      amp::unordered_set<int, Hasher::IntHashFcn, IntCompare> set;

      //for (int k = 0; k < 23; ++k) {
      //   auto b = set.GetBucket(k);
      //   printf("%d b");
      //}

      int maxNum = 100;
      set.insert(maxNum);
      for (int k = 0; k < maxNum; ++k) {
         set.insert(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set.Exists(k)) {
            printf("number not found: k\n", k);
         }
      }

      int num1 = 50, num2 = 100;

      set.erase(num1);
      set.erase(num2);

      if (set.Exists(num1)) {
         printf("not removing elements properly: %d\n", num1);
      }
      if (set.Exists(num2)) {
         printf("not removing elements properly: %d\n", num2);
      }

      for (int k = 0; k < maxNum - 20; ++k) {
         set.erase(k);
         if (set.Exists(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      set.clear();

      for (int k = 0; k < maxNum; ++k) {
         if (set.Exists(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      printf("SetTest::TestInt() completed.\n");
   }

   __device__ void SetTest::TestInt2()
   {
      amp::unordered_set<int, Hasher::IntHashFcn, IntCompare> set1;

      int maxNum = 100;

      for (int k = 0; k < maxNum; ++k) {
         set1.insert(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set1.Exists(k)) {
            printf("number not found: k\n", k);
         }
      }

      amp::unordered_set<int, Hasher::IntHashFcn, IntCompare> set2;
      for (int k = 0; k < maxNum; ++k) {
         set2.insert(k);
      }

      unsigned int h1 = set1.HashCode();
      unsigned int h2 = set2.HashCode();
      if (h1 != h2) {
         printf("HashCodes are different: %u, %u", h1, h2);
      }

      printf("SetTest::TestInt2() completed.\n");
   }

   typedef amp::unordered_set<String::String, String::HashFcn, String::CompareFcn> StringTestT;
   typedef amp::Iterator<String::String, String::HashFcn, String::CompareFcn> StringTestTItr;

   __device__ void SetTest::TestString()
   {
      StringTestT set;
      
      String::String a("a");
      String::String b("b");
      String::String a2("a");
      String::String a3("a");
      String::String abc("abc");

      String::CompareFcn cmp;
      if (cmp(a, b)) {
         printf("wrong: %s, %s\n", a.Get(), b.Get());
      }

      set.insert(a);
      set.insert(b);
      set.insert(abc);

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

      StringTestTItr itr(set);
      while (itr.HasNext()) {
         //printf("%s ", itr.GetValue().Get());
         itr.Next();
      }
      //printf("\n");

      StringTestT set2 = set;
      String::String xyz("XYZ");
      set2.insert(xyz);

      if (!(set == set2)) {
         printf("sets are different: %d, %d\n", set.HashCode(), set2.HashCode());
      }

      int c2 = 0;
      StringTestTItr itr2(set2);
      while (itr2.HasNext()) {
         //printf("%s ", itr2.GetValue().Get());
         ++c2;
         itr2.Next();
      }
      //printf("\n");
      if (c2 != set2.size()) {
         printf("different sizes: (%d, %d)", c2, set2.size());
      }

      printf("SetTest::TestString() completed.\n");
   }

   //class ClassSet
   //{
   //public:
   //   __host__ __device__ ClassSet() {}
   //   __host__ __device__ ClassSet(StringTestT& s) : s(s) {}
   //   __host__ __device__ StringTestT& GetSet() { return s; }

   //private:
   //   StringTestT s;
   //};

   //struct ClassSetCompare
   //{
   //   __host__ __device__ inline bool operator() (ClassSet& lhs, ClassSet& rhs)
   //   {
   //      return lhs.GetSet() == rhs.GetSet();
   //   }
   //};

   //struct ClassSetHash
   //{
   //   __host__ __device__ inline unsigned int operator() (ClassSet& lhs)
   //   {
   //      return lhs.GetSet().HashCode();
   //   }
   //};

   //typedef Set::Set<ClassSet, ClassSetCompare, ClassSetHash> SetOfSetT;
   //typedef Set::Iterator<ClassSet, ClassSetCompare, ClassSetHash> SetOfSetTItr;

   //__device__ void SetTest::TestSetOfSetOfString()
   //{
   //   StringTestT set0, set1;
   //   String::String str0("A"), str1("B"), str2("C");
   //   set0.insert(str0);
   //   set0.insert(str1);
   //   set1.insert(str2);

   //   ClassSet cs0(set0);
   //   ClassSet cs1(set1);

   //   SetOfSetT scs;
   //   scs.insert(cs0);
   //   scs.insert(cs1);
   //   //printf("scs: %d\n", scs.Size());

   //   SetOfSetTItr itr1(scs);
   //   while (itr1.HasNext()) {
   //      ClassSet cstmp = itr1.GetValue();
   //      //printf("itr1: %d\n", cstmp.GetSet().Size());

   //      StringTestT tmpset = cstmp.GetSet();
   //      StringTestTItr itr2(tmpset);
   //      while (itr2.HasNext()) {
   //         //printf("%s ", itr2.GetValue().Get());
   //         itr2.Next();
   //      }
   //      //printf("\n");
   //      itr1.Next();
   //   }

   //   StringTestT set2(set0);
   //   ClassSet cs2(set2);
   //   StringTestT tmpset = cs2.GetSet();
   //   printf("cs2: %d\n", tmpset.Size());
   //   StringTestTItr itr(tmpset);
   //   while (itr.HasNext()) {
   //      printf("%s ", itr.GetValue().Get());
   //      itr.Next();
   //   }
   //   printf("\n");

   //   if (!scs.Exists(cs2)) {
   //      printf("failed to check for existence by value\n");
   //   }

   //   if (scs.HashCode() != 2147483668) {
   //      printf("hash codes not equal: %u %u\n", scs.HashCode(), 2147483668);
   //   }

   //   printf("SetTest::TestSetOfSetOfString() completed.\n");
   //}
}
