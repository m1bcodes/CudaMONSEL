#include "Amphibian/Tests/LinkedListTest.cuh"
#include "Amphibian/string.cuh"

#include <stdio.h>

namespace LinkedListTest
{
   class TestClassA
   {
   public:
      __host__ __device__ TestClassA() : val(0), list(nullptr)
      {
      }

      __host__ __device__ TestClassA(double v, char name[], double s) : val(v)
      {
         amp::string nameStr(name);
         LinkedListKV::InsertHead<amp::string, double>(&list, nameStr, s);
      }

   private:
      double val;
      LinkedListKV::Node<amp::string, double>* list;
   };

   __host__ __device__ static void assertTrue(bool expr)
   {
      if (!expr) printf("wrong!\n");
   }

   __host__ __device__ void PrintList(LinkedList::Node<LinkedListTestType>* head)
   {
      while (true) {
         if (head != nullptr) {
            printf("%d, ", head->GetValue());
         }
         else {
            printf("\n");
            break;
         }
         head = head->GetNext();
      }
   }

   __host__ __device__ LinkedListTest::LinkedListTest()
   {
   }

   __host__ __device__ void LinkedListTest::InsertionTest()
   {
      LinkedList::Node<int> * head = nullptr;
      LinkedList::Node<int>** newHeadAddr = &head;
      int v1 = 0, v2 = 1, v3 = 2;
      LinkedList::InsertNext(newHeadAddr, v1);
      LinkedList::InsertNext(newHeadAddr, v2);
      newHeadAddr = (*newHeadAddr)->GetNextAddr();
      LinkedList::InsertNext(newHeadAddr, v3);
   
      PrintList(head);
      LinkedList::Node<int> * head2 = nullptr;
      LinkedList::DeepCopy(&head2, head);
      LinkedList::RemoveAll(&head);
      PrintList(head);
      PrintList(head2);
      LinkedList::RemoveAll(&head2);
      PrintList(head2);
   }

   __host__ __device__ void LinkedListTest::TestAddAllAsSet()
   {
      LinkedList::Node<LinkedListTestType>* head = nullptr;
      int a[] = { 0, 1, 2, 3 };
      int b[] = { 0, 1, 2, 3 };
      int c[] = { 5 };

      LinkedList::Node<LinkedListTestType>* alist = nullptr, * blist = nullptr, * clist = nullptr;
      LinkedList::BuildList(&alist, a, 4);
      LinkedList::BuildList(&blist, b, 4);
      LinkedList::BuildList(&clist, c, 1);
      printf("LinkedListTest::TestAddAllAsSet completed\n");
   }

   __host__ __device__ void BuildList1(LinkedListKV::Node<amp::string, float>** head)
   {
      amp::string a("a");
      amp::string b("b");
      amp::string c("c");
      float v1 = 0.0f;
      float v2 = 1.0f;
      float v3 = 2.0f;
      LinkedListKV::InsertHead(head, a, v1);
      LinkedListKV::InsertHead(head, b, v2);
      LinkedListKV::InsertHead(head, c, v3);
   }
   
   __host__ __device__ void BuildList2(LinkedListKV::Node<amp::string, float>** head)
   {
      amp::string a("a");
      amp::string b("b");
      amp::string a1("a");
      float v1 = 0.0f;
      float v2 = 1.0f;
      float v3 = .1f;
      LinkedListKV::InsertHead(head, a, v1);
      LinkedListKV::InsertHead(head, b, v2);
      LinkedListKV::InsertHead(head, a1, v3);
   }
   
   __host__ __device__ void PrintListInOrder(LinkedListKV::Node<amp::string, float>* head)
   {
      if (head != nullptr) {
         printf("%s: %f\n", head->GetKey().c_str(), head->GetValue());
         PrintListInOrder(head->GetNext());
      }
   }
   
   __global__ void SetTest(amp::pStrCmp h_strCmpPointFunction)
   {
   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      __syncthreads();
   #endif

      LinkedListKV::Node<amp::string, float>* head1 = nullptr;
      LinkedListKV::Node<amp::string, float>* head2 = nullptr;
      printf("A\n");
      BuildList1(&head1);
      BuildList2(&head2);
      printf("B\n");
      PrintListInOrder(head1);
      PrintListInOrder(head2);
   
   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      __syncthreads();
   #endif
   }

   //__host__ __device__ amp::pAreEqual pEqual = amp::AreEqual;
   //
   //int main()
   //{
   //   amp::pAreEqual h_1;
   //   cudaMemcpyFromSymbol(&h_1, pEqual, sizeof(amp::pAreEqual));
   //
   //   kernel << <1, 1 >> >(h_1);
   //   checkCudaErrors(cudaDeviceSynchronize());
   //   checkCudaErrors(cudaGetLastError());
   //
   //   return 0;
   //}

   __host__ __device__ void TestListKV()
   {
      LinkedListKV::Node<amp::string, int>* pairlist = new LinkedListKV::Node<amp::string, int>("a", 1, nullptr);
      amp::string k = "b";
      int v = 2;
      LinkedListKV::InsertHead(&pairlist, k, v);
      k = "c";
      v = 3;
      LinkedListKV::InsertHead(&pairlist, k, v);

      while (pairlist) {
         printf("%s, %d\n", ((amp::string&)(pairlist->first)).c_str(), (int)pairlist->second);
         pairlist = pairlist->GetNext();
      }

      printf("LinkedListTest::TestListKV() completed.\n");
   }

   __host__ __device__ void testDLinkedList()
   {
      DLinkedList::Node<int>* head = nullptr;
      int a = 0;
      DLinkedList::Node<int>* newone0 = new DLinkedList::Node<int>(nullptr, nullptr, a);
      DLinkedList::Node<int>::Insert(&head, newone0);
      a = 1;
      DLinkedList::Node<int>* newone1 = new DLinkedList::Node<int>(nullptr, nullptr, a);
      DLinkedList::Node<int>::Insert(&head, newone1);
      a = 2;
      DLinkedList::Node<int>* newone2 = new DLinkedList::Node<int>(nullptr, nullptr, a);
      DLinkedList::Node<int>::Insert(&head, newone2);

      auto itr0 = newone0;
      while (itr0) {
         printf("%d ", itr0->GetValue());
         itr0 = itr0->GetPrev();
      }
      printf("\n");

      auto itr1 = newone2;
      while (itr1) {
         printf("%d ", itr1->GetValue());
         itr1 = itr1->GetNext();
      }
      printf("\n");

      DLinkedList::Node<int>::Remove(&head);
      DLinkedList::Node<int>::Remove(&head);
      DLinkedList::Node<int>::Remove(&head);
      DLinkedList::Node<int>::Remove(&head);
      assertTrue(head == nullptr);

      //int i;
      //for (i = 0; i < 100000; ++i) {
      //   DLinkedList::Node<int>* newone = new DLinkedList::Node<int>(nullptr, nullptr, i);
      //   DLinkedList::Node<int>::Insert(&head, newone);
      //}
      //assertTrue(i == 100000);

      //for (i = 0; i < 100000; ++i) {
      //   if (head == nullptr) break;
      //   DLinkedList::Node<int>::Remove(&head);
      //}
      //assertTrue(head == nullptr);
      //assertTrue(i == 100000);

      printf("LinkedListTest::testDLinkedList() completed.\n");
   }
}