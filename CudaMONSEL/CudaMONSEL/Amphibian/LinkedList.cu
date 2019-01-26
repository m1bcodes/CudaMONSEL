#ifndef LINKED_LIST_CU
#define LINKED_LIST_CU

#include "cuda_runtime.h"

template<typename Key, typename Value>
class Node
{
public:
   __host__ __device__ static void InsertHead(Node** head, Key key, Value val);
   __host__ __device__ static void RemoveHead(Node** head);

   __host__ __device__ static void RemoveAll(Node** head);

   __host__ __device__ Key GetKey();
   __host__ __device__ Value GetValue();
   __host__ __device__ Node<Key, Value>* GetNext();

private:
   Key key;
   Value val;
   Node* next;
};

#endif