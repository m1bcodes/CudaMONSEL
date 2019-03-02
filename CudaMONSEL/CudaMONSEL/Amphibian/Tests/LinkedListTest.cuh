#ifndef _LINKED_LIST_TEST_CUH_
#define _LINKED_LIST_TEST_CUH_

#include "..\LinkedList.cuh"

#include <cuda_runtime.h>

namespace LinkedListTest
{
   typedef int LinkedListTestType;

   __device__ void PrintList(LinkedList::Node<LinkedListTestType>* head);

   class LinkedListTest
   {
   public:
      __device__ LinkedListTest();
      __device__ void InsertionTest();
      __device__ void TestAddAllAsSet();
   };
}

#endif