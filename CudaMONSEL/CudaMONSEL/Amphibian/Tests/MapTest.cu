//#include "Amphibian/Tests/MapTest.cuh"
//
//#include "Amphibian/String.cuh"
//
//#include <math.h>
//
//namespace MapTest
//{
//   __device__ void AssertEqual(int a, int b)
//   {
//      if (a != b) {
//         printf("not equal: (%d, %d)\n", a, b);
//      }
//   }
//
//   __device__ MapTest::MapTest()
//   {
//   }
//
//   typedef Map::Map<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn> IntTestType;
//   typedef Map::Iterator<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn> IntTestTypeItr;
//
//   __device__ void MapTest::TestInteger()
//   {
//      int k = 0, v = 1;
//      IntTestType m1 =
//         CreateMapA<int, int, Comparator::IntCompareFcn, Comparator::IntCompareFcn, Hasher::IntHashFcn, Hasher::IntHashFcn>(k, v);
//      int k0 = 1, v0 = 1;
//      m1.Put(k0, v0);
//
//      AssertEqual(m1.Size(), 2);
//
//      IntTestType m2(m1);
//      AssertEqual(m2.Size(), 2);
//      int k1 = 2, v1 = 3;
//      m2.Put(k1, v1);
//      AssertEqual(m2.Size(), 3);
//      AssertEqual(m1.Size(), 2);
//
//      int c1 = 0;
//      IntTestTypeItr itr1(m1);
//      while (itr1.HasNext()) {
//         //printf("(%d, %d) ", itr1.GetKey(), itr1.GetValue());
//         c1++;
//         itr1.Next();
//      }
//      AssertEqual(c1, 2);
//
//      int c2 = 0;
//      IntTestTypeItr itr2(m2);
//      while (itr2.HasNext()) {
//         //printf("(%d, %d) ", itr2.GetKey(), itr2.GetValue());
//         c2++;
//         itr2.Next();
//      }
//      AssertEqual(c2, 3);
//
//      IntTestType m3 = m2;
//      AssertEqual(m3.Size(), 3);
//      int k2 = 2, v2 = 3;
//      m3.Put(k2, v2);
//      AssertEqual(m3.Size(), 3);
//      unsigned int h2 = m2.HashCode();
//      unsigned int h3 = m3.HashCode();
//      if (h2 != h3) {
//         printf("HashCodes are different: %d, %d\n", h2, h3);
//      }
//
//      printf("MapTest::TestInteger() completed\n");
//   }
//
//   typedef Map::Map<String::String, double, String::CompareFcn, Comparator::DoubleCompareFcn, String::HashFcn, Hasher::DoubleHashFcn> StringTestType;
//   typedef Map::Iterator<String::String, double, String::CompareFcn, Comparator::DoubleCompareFcn, String::HashFcn, Hasher::DoubleHashFcn> StringTestTypeItr;
//
//   __device__ StringTestType makeMap()
//   {
//      StringTestType map;
//      double v0 = 1, v1 = 2, v2 = 3, v3 = 4, v4 = 5;
//      String::String A("A"), B("B"), C("C"), D("D"), E("E");
//      map.Put(A, v0);
//      map.Put(B, v1);
//      map.Put(C, v2);
//      map.Put(D, v3);
//      map.Put(E, v4);
//      return map;
//   }
//
//   __device__ void MapTest::TestString()
//   {
//      StringTestType m1;
//      double v0 = 1, v1 = 2, v2 = 3, v3 = 4;
//      String::String A("A"), B("B"), C("C"), D("D"), E("E");
//      m1.Put(A, v0);
//      m1.Put(B, v1);
//
//      AssertEqual(m1.Size(), 2);
//
//      StringTestType m2(m1);
//      AssertEqual(m2.Size(), 2);
//      m2.Put(C, v2);
//      AssertEqual(m2.Size(), 3);
//      AssertEqual(m1.Size(), 2);
//
//      int c1 = 0;
//      StringTestTypeItr itr1(m1);
//      while (itr1.HasNext()) {
//         //printf("(%s, %lf) ", itr1.GetKey().Get(), itr1.GetValue());
//         c1++;
//         itr1.Next();
//      }
//      AssertEqual(c1, 2);
//      //printf("\n");
//
//      int c2 = 0;
//      StringTestTypeItr itr2(m2);
//      while (itr2.HasNext()) {
//         //printf("(%s, %lf) ", itr2.GetKey().Get(), itr2.GetValue());
//         c2++;
//         itr2.Next();
//      }
//      AssertEqual(c2, 3);
//      //printf("\n");
//
//      StringTestType m3 = m2;
//      AssertEqual(m3.Size(), 3);
//      m3.Put(D, v3);
//      AssertEqual(m3.Size(), 4);
//
//      int c3 = 0;
//      StringTestTypeItr itr3(m3);
//      while (itr3.HasNext()) {
//         //printf("(%s, %lf) ", itr3.GetKey().Get(), itr3.GetValue());
//         c3++;
//         itr3.Next();
//      }
//      AssertEqual(m3.Size(), 4);
//      //printf("\n");
//
//      StringTestType map4;
//      double v5 = ::sqrt(0.05), v6 = ::sqrt(0.04);
//      String::String V1("V1"), V2("V2");
//      map4.Put(V1, v5);
//      map4.Put(V2, v6);
//
//      int c4 = 0;
//      StringTestTypeItr itr4(map4);
//      while (itr4.HasNext()) {
//         //printf("(%s, %lf) ", itr4.GetKey().Get(), itr4.GetValue());
//         c4++;
//         itr4.Next();
//      }
//      AssertEqual(map4.Size(), 2);
//      AssertEqual(c4, 2);
//
//      StringTestType m5 = map4;
//      if (!(m5 == map4)) {
//         printf("maps are different\n");
//      }
//      if (!(m5.Size() == map4.Size())) {
//         printf("maps sizes are different: %d, %d\n", map4.Size(), m5.Size());
//      }
//      if (m5.HashCode() != map4.HashCode()) {
//         printf("maps hashcodes are different\n");
//      }
//
//      StringTestTypeItr itr5(m5);
//      while (itr5.HasNext()) {
//         //printf("(%s, %lf) ", itr5.GetKey().Get(), itr5.GetValue());
//         itr5.Next();
//      }
//
//      StringTestType m6 = makeMap();
//      StringTestTypeItr itr6(m6);
//      while (itr6.HasNext()) {
//         //printf("(%s, %lf) ", itr6.GetKey().Get(), itr6.GetValue());
//         itr6.Next();
//      }
//      //printf("\n");
//
//      printf("MapTest::TestString() completed\n");
//   }
//
//   class TestClass
//   {
//   public:
//      __device__ TestClass() {}
//      __device__ StringTestType& GetMap() { return m; }
//
//   private:
//      StringTestType m;
//   };
//
//   struct TestClassCompare
//   {
//      __device__ inline bool operator() (TestClass& lhs, TestClass& rhs)
//      {
//         return lhs.GetMap() == rhs.GetMap();
//      }
//   };
//
//   struct TestClassHashFcn
//   {
//      __device__ inline unsigned int operator() (TestClass& t)
//      {
//         return t.GetMap().HashCode();
//      }
//   };
//
//   typedef Map::Map<String::String, TestClass, String::CompareFcn, TestClassCompare, String::HashFcn, TestClassHashFcn> MapTestT;
//   typedef Map::Iterator<String::String, TestClass, String::CompareFcn, TestClassCompare, String::HashFcn, TestClassHashFcn> MapTestTItr;
//
//   __device__ void MapTest::TestMapOfMap()
//   {
//      MapTestT m1;
//      TestClass a;
//      TestClass b;
//
//      double v1 = 1.2, v2 = 1.2, v3 = 3.3;
//      String::String s0("A1"), s1("A2"), s2("A3"), s3("B1");
//      String::String s4("A"), s5("B");
//      a.GetMap().Put(s0, v1);
//      a.GetMap().Put(s1, v2);
//      a.GetMap().Put(s2, v3);
//      b.GetMap().Put(s3, v1);
//      m1.Put(s4, a);
//      m1.Put(s5, b);
//
//      printf("MapTest::TestMapOfMap() completed\n");
//   }
//
//   __device__ void MapTest::TestAggregate()
//   {
//      StringTestType t;
//      double v0 = 0.1, v1 = 0.5, v2 = 100;
//      String::String s0("A"), s1("B"), s2("C");
//      t.Put(s0, v0);
//      t.Put(s1, v1);
//      t.Put(s2, v2);
//
//      double ret1 = t.Aggregate([](double a) {return a*a; });
//      double ret2 = t.Aggregate([](double a) {return a*a; });
//      double ret3 = t.Aggregate([](double a) {return a*a; });
//      double ans = 0.1*0.1 + 0.5*0.5 + 100 * 100;
//      if (ret1 != ans) {
//         printf("wrong : %lf, %lf", ret1, ans);
//      }
//      if (ret2 != ans) {
//         printf("wrong : %lf, %lf", ret2, ans);
//      }
//      if (ret3 != ans) {
//         printf("wrong : %lf, %lf", ret3, ans);
//      }
//      printf("MapTest::TestAggregate() completed\n");
//   }
//}
