#include "Amphibian/Tests/LinkedListTest.cuh"
#include "Amphibian/string.cuh"

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
         amp::string nameStr(name);
         LinkedListKV::InsertHead<amp::string, double>(&list, nameStr, s);
      }

      double val;
      LinkedListKV::Node<amp::string, double>* list;
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

   __device__ void BuildList1(LinkedListKV::Node<amp::string, float>** head)
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
   
   __device__ void BuildList2(LinkedListKV::Node<amp::string, float>** head)
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
   
   __device__ void PrintListInOrder(LinkedListKV::Node<amp::string, float>* head)
   {
      if (head != NULL) {
         printf("%s: %f\n", head->GetKey().c_str(), head->GetValue());
         PrintListInOrder(head->GetNext());
      }
   }
   
   __global__ void SetTest(amp::pStrCmp h_strCmpPointFunction)
   {
   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      __syncthreads();
   #endif

      LinkedListKV::Node<amp::string, float>* head1 = NULL;
      LinkedListKV::Node<amp::string, float>* head2 = NULL;
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

   //__device__ amp::pAreEqual pEqual = amp::AreEqual;
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
}