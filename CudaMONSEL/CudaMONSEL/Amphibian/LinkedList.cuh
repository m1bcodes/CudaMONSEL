#ifndef LINKED_LIST_CUH
#define LINKED_LIST_CUH

#include "cuda_runtime.h"

template<typename Key, typename Value>
class Node
{
public:
   __host__ __device__ static void InsertHead(Node** head, Key key, Value val);
   __host__ __device__ static void RemoveHead(Node** head);

   __host__ __device__ static void RemoveAll(Node** head);

   __host__ __device__ Node();
   __host__ __device__ Node(Key, Value, Node *);

   __host__ __device__ Key GetKey();
   __host__ __device__ Value GetValue();
   __host__ __device__ Node<Key, Value>* GetNext();

private:
   Key key;
   Value val;
   Node* next;
};


template<typename Key, typename Value>
__host__ __device__ void Node<Key, Value>::InsertHead(Node<Key, Value>** head, Key k, Value v)
{
   Node<Key, Value>* newOne = new Node<Key, Value>(k, v, NULL);

   if (*head == NULL) {
      *head = newOne;
   }
   else {
      newOne->next = *head;
      *head = newOne;
   }
}

template<typename Key, typename Value>
__host__ __device__ void Node<Key, Value>::RemoveHead(Node<Key, Value>** head)
{
   if (head == NULL) {
      return;
   }

   Node<Key, Value>* tmp = (*head)->next;
   delete *head;
   *head = tmp;
}

template<typename Key, typename Value>
__host__ __device__ void Node<Key, Value>::RemoveAll(Node<Key, Value>** head)
{
   while (*head != NULL) {
      Node<Key, Value>::RemoveHead(head);
   }
}

template<typename Key, typename Value>
__host__ __device__ Node<Key, Value>::Node<Key, Value>()
{
}

template<typename Key, typename Value>
__host__ __device__ Node<Key, Value>::Node<Key, Value>(Key k, Value v, Node* n) : key(k), val(v), next(n)
{
}

template<typename Key, typename Value>
__host__ __device__ Key Node<Key, Value>::GetKey()
{
   return key;
}

template<typename Key, typename Value>
__host__ __device__ Value Node<Key, Value>::GetValue()
{
   return val;
}

template<typename Key, typename Value>
__host__ __device__ Node<Key, Value>* Node<Key, Value>::GetNext()
{
   return next;
}

#endif