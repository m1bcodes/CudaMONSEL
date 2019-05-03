#include "LinkedListTest.cuh"
#include "Amphibian\String.cuh"

#include <stdio.h>

namespace LinkedListTest
{
   class TestClassA
   {
   public:
      __device__ TestClassA() : val(0), list(NULL)
      {
      }

      __device__ TestClassA(double v, char name[], double s) : val(v)
      {
         String::String nameStr(name);
         LinkedListKV::InsertHead<String::String, double>(&list, nameStr, s);
      }

      double val;
      LinkedListKV::Node<String::String, double>* list;
   };

   __device__ void PrintList(LinkedList::Node<LinkedListTestType>* head)
   {
      while (true) {
         if (head != NULL) {
            printf("%d, ", head->GetValue());
         }
         else {
            printf("\n");
            break;
         }
         head = head->GetNext();
      }
   }

   __device__ LinkedListTest::LinkedListTest()
   {
   }

   __device__ void LinkedListTest::InsertionTest()
   {
      LinkedList::Node<int> * head = NULL;
      LinkedList::Node<int>** newHeadAddr = &head;
      int v1 = 0, v2 = 1, v3 = 2;
      LinkedList::InsertNext(newHeadAddr, v1);
      LinkedList::InsertNext(newHeadAddr, v2);
      newHeadAddr = (*newHeadAddr)->GetNextAddr();
      LinkedList::InsertNext(newHeadAddr, v3);
   
      PrintList(head);
      LinkedList::Node<int> * head2 = NULL;
      LinkedList::DeepCopy(&head2, head);
      LinkedList::RemoveAll(&head);
      PrintList(head);
      PrintList(head2);
      LinkedList::RemoveAll(&head2);
      PrintList(head2);
   }

   __device__ void LinkedListTest::TestAddAllAsSet()
   {
      LinkedList::Node<LinkedListTestType>* head = NULL;
      int a[] = { 0, 1, 2, 3 };
      int b[] = { 0, 1, 2, 3 };
      int c[] = { 5 };

      LinkedList::Node<LinkedListTestType>* alist = NULL, * blist = NULL, * clist = NULL;
      LinkedList::BuildList(&alist, a, 4);
      LinkedList::BuildList(&blist, b, 4);
      LinkedList::BuildList(&clist, c, 1);
      printf("LinkedListTest::testAddAllAsSet completed\n");
   }

   __device__ void BuildList1(LinkedListKV::Node<String::String, float>** head)
   {
      String::String a("a");
      String::String b("b");
      String::String c("c");
      float v1 = 0.0f;
      float v2 = 1.0f;
      float v3 = 2.0f;
      LinkedListKV::InsertHead(head, a, v1);
      LinkedListKV::InsertHead(head, b, v2);
      LinkedListKV::InsertHead(head, c, v3);
   }
   
   __device__ void BuildList2(LinkedListKV::Node<String::String, float>** head)
   {
      String::String a("a");
      String::String b("b");
      String::String a1("a");
      float v1 = 0.0f;
      float v2 = 1.0f;
      float v3 = .1f;
      LinkedListKV::InsertHead(head, a, v1);
      LinkedListKV::InsertHead(head, b, v2);
      LinkedListKV::InsertHead(head, a1, v3);
   }
   
   __device__ void PrintListInOrder(LinkedListKV::Node<String::String, float>* head)
   {
      if (head != NULL) {
         printf("%s: %f\n", head->GetKey().Get(), head->GetValue());
         PrintListInOrder(head->GetNext());
      }
   }
   
   __global__ void SetTest(String::pStrCmp h_strCmpPointFunction)
   {
   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      __syncthreads();
   #endif

      //LinkedListKV::Node<String::String, float>* head1 = NULL;
      //LinkedListKV::Node<String::String, float>* head2 = NULL;
      //printf("A\n");
      //BuildList1(&head1);
      //BuildList2(&head2);
      //printf("B\n");
      //PrintListInOrder(head1);
      //PrintListInOrder(head2);
      //printf("C\n");
      //printf("%d\n", LinkedListKV::IsSet<String::String, float>(head1, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
      //printf("%d\n", LinkedListKV::IsSet<String::String, float>(head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
      //printf("D\n");
      //printf("%d\n", LinkedListKV::AreEquivalentSets<String::String, float>(head1, head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
      //LinkedListKV::RemoveRepeatedNodes<String::String, float>(&head2, h_strCmpPointFunction, [](float a, float b) { return a == b; });
      //PrintListInOrder(head2);
      //printf("E\n");
      //LinkedListKV::Remove<String::String, float>(&head2, String::String("a"), h_strCmpPointFunction);
      //PrintListInOrder(head2);
      //printf("%d\n", LinkedListKV::AreEquivalentSets<String::String, float>(head1, head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
      //LinkedListKV::RemoveAll(&head1);
      //LinkedListKV::RemoveAll(&head2);
   
   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      __syncthreads();
   #endif
   }

   //__device__ String::pAreEqual pEqual = String::AreEqual;
   //
   //int main()
   //{
   //   String::pAreEqual h_1;
   //   cudaMemcpyFromSymbol(&h_1, pEqual, sizeof(String::pAreEqual));
   //
   //   kernel << <1, 1 >> >(h_1);
   //   checkCudaErrors(cudaDeviceSynchronize());
   //   checkCudaErrors(cudaGetLastError());
   //
   //   return 0;
   //}
}