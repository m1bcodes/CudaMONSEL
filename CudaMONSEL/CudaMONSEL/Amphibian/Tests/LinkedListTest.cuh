#ifndef _LINKED_LIST_TEST_CUH_
#define _LINKED_LIST_TEST_CUH_

#include "Amphibian\LinkedList.cuh"

#include <cuda_runtime.h>

namespace LinkedListTest
{
   typedef int LinkedListTestType;

   __host__ __device__ void PrintList(LinkedList::Node<LinkedListTestType>* head);

   class LinkedListTest
   {
   public:
      __host__ __device__ LinkedListTest();
      __host__ __device__ void InsertionTest();
      __host__ __device__ void TestAddAllAsSet();
   };

   __host__ __device__ void TestListKV();

   __host__ __device__ void testDLinkedList();
}

#endif