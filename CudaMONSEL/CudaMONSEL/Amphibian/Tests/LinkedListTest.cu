#include "LinkedListTest.cuh"

#include <stdio.h>

namespace LinkedListTest
{
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

      LinkedList::AddAllAsSet<LinkedListTestType>(&head, alist, [](LinkedListTestType a, LinkedListTestType b) { return a == b; });
      PrintList(head);
      LinkedList::AddAllAsSet<LinkedListTestType>(&head, blist, [](LinkedListTestType a, LinkedListTestType b) { return a == b; });
      PrintList(head);
      LinkedList::AddAllAsSet<LinkedListTestType>(&head, clist, [](LinkedListTestType a, LinkedListTestType b) { return a == b; });
      PrintList(head);
      printf("LinkedListTest::testAddAllAsSet completed\n");
   }
}