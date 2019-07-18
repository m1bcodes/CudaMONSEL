#include "Amphibian/Tests/SetTest.cuh"

#include "Amphibian/unordered_set.cuh"
#include "Amphibian/String.cuh"
#include "Amphibian/Hasher.cuh"

#include <stdio.h>

namespace SetTest
{
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

   __host__ __device__ SetTest::SetTest()
   {
   }

   __host__ __device__ void SetTest::testIntBasic()
   {
      Comparator::IntCompareFcn cmp;
      int a = 1, b = 1;
      if (!cmp(a, b)) {
         printf("bad: %d, %d\n", a, b);
      }
      amp::unordered_set<int, Hasher::IntHashFcn, Comparator::IntCompareFcn> set;
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

      printf("SetTest::testIntBasic() completed.\n");
   }

   __host__ __device__ void SetTest::testInt()
   {
      amp::unordered_set<int, Hasher::IntHashFcn, Comparator::IntCompareFcn> set;

      //for (int k = 0; k < 23; ++k) {
      //   auto b = set.GetBucket(k);
      //   printf("%d b");
      //}

      const int maxNum = 100;
      set.insert(maxNum);
      for (int k = 0; k < maxNum; ++k) {
         set.insert(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set.contains(k)) {
            printf("number not found: k\n", k);
         }
      }

      const int num1 = 50, num2 = 100;

      set.erase(num1);
      set.erase(num2);

      if (set.contains(num1)) {
         printf("not removing elements properly: %d\n", num1);
      }
      if (set.contains(num2)) {
         printf("not removing elements properly: %d\n", num2);
      }

      for (int k = 0; k < maxNum - 20; ++k) {
         set.erase(k);
         if (set.contains(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      set.clear();

      for (int k = 0; k < maxNum; ++k) {
         if (set.contains(k)) {
            printf("not removing number: %d\n", k);
         }
      }

      printf("SetTest::testInt() completed.\n");
   }

   __host__ __device__ void SetTest::testInt2()
   {
      amp::unordered_set<int, Hasher::IntHashFcn, Comparator::IntCompareFcn> set1;

      int maxNum = 100;

      for (int k = 0; k < maxNum; ++k) {
         set1.insert(k);
      }
      for (int k = 0; k < maxNum; ++k) {
         if (!set1.contains(k)) {
            printf("number not found: k\n", k);
         }
      }
      for (auto k : set1) {
         printf("%d ", k);
      }
      printf("\n");

      amp::unordered_set<int, Hasher::IntHashFcn, Comparator::IntCompareFcn> set2;
      for (int k = 0; k < maxNum; ++k) {
         set2.insert(k);
      }

      unsigned int h1 = set1.hashCode();
      unsigned int h2 = set2.hashCode();
      if (h1 != h2) {
         printf("HashCodes are different: %u, %u", h1, h2);
      }

      printf("SetTest::testInt2() completed.\n");
   }

   typedef amp::unordered_set<amp::string, amp::string_hash, amp::string_cmp> StringTestT;
   typedef amp::unordered_set<amp::string, amp::string_hash, amp::string_cmp>::const_iterator StringTestTItr;

   __host__ __device__ void SetTest::testString()
   {
      StringTestT set;
      
      amp::string a("a");
      amp::string b("b");
      amp::string a2("a");
      amp::string a3("a");
      amp::string abc("abc");

      amp::string_cmp cmp;
      if (cmp(a, b)) {
         printf("wrong: %s, %s\n", a.c_str(), b.c_str());
      }

      set.insert(a);
      set.insert(b);
      set.insert(abc);

      if (!set.contains(a2)) {
         printf("does not exist: %s\n", a2.c_str());
      }

      if (!set.contains(a3)) {
         printf("does not exist: %s\n", a2.c_str());
      }

      amp::string d("d");

      if (set.contains(d)) {
         printf("exists: %s\n", d.c_str());
      }

      StringTestTItr itr(set);
      while (itr.HasNext()) {
         //printf("%s ", itr.GetValue().c_str());
         itr.next();
      }
      //printf("\n");

      for (auto str : set) {
         printf("%s ", str.c_str());
      }
      printf("\n");

      StringTestT set2 = set;
      amp::string xyz("XYZ");
      set2.insert(xyz);

      if (!(set == set2)) {
         printf("sets are different: %d, %d\n", set.hashCode(), set2.hashCode());
      }

      int c2 = 0;
      StringTestTItr itr2(set2);
      while (itr2.HasNext()) {
         //printf("%s ", itr2.GetValue().c_str());
         ++c2;
         itr2.next();
      }
      //printf("\n");
      if (c2 != set2.size()) {
         printf("different sizes: (%d, %d)", c2, set2.size());
      }

      printf("SetTest::testString() completed.\n");
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
   //   amp::string str0("A"), str1("B"), str2("C");
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
   //         //printf("%s ", itr2.GetValue().c_str());
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
   //      printf("%s ", itr.GetValue().c_str());
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
