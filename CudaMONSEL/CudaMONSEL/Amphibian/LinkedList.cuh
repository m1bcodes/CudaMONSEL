#ifndef LINKED_LIST_CUH
#define LINKED_LIST_CUH

#include "cuda_runtime.h"

template<typename Key, typename Value>
class Node
{
public:
   __host__ __device__ static void InsertHead(Node** head, Key key, Value val);
   __host__ __device__ static void RemoveHead(Node** head);

   __host__ __device__ static void Remove(Node** head, Key k, bool equals(Key, Key));
   __host__ __device__ static void RemoveAll(Node** head);

   __host__ __device__ static bool AreEquivalentNodes(Node*, Node*, bool equalKeys(Key, Key), bool equalValues(Value, Value));
   __host__ __device__ static bool IsSet(Node*, bool equalKeys(Key, Key), bool equalValues(Value, Value));
   __host__ __device__ static bool AreEquivalentSets(Node*, Node*, bool equalKeys(Key, Key), bool equalValues(Value, Value));

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
__host__ __device__ static void Remove(Node<Key, Value>** head, Key k, bool equals(Key, Key))
{
   while (*head != NULL) {
      if (equals((*head)->GetKey(), k)) {
         Node<Key, Value>::RemoveHead(head);
      }
      else {
         Node<Key, Value>::Remove((*head)->GetNext(), k, equals);
      }
   }
}

template<typename Key, typename Value>
__host__ __device__ static bool Node<Key, Value>::AreEquivalentNodes(Node<Key, Value>* head1, Node<Key, Value>* head2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   if ((*head1) == NULL && (*head2) == NULL) {
      return true;
   }

   if ((*head1) == NULL || (*head2) == NULL) {
      return false;
   }

   return (equalKeys(head1->GetKey(), head2->GetKey()) && equalValues(head1->GetValue(), head2->GetValue()));
}

template<typename Key, typename Value>
__host__ __device__ static bool Node<Key, Value>::IsSet(Node<Key, Value>* head, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   return false;
}

template<typename Key, typename Value>
__host__ __device__ static bool Node<Key, Value>::AreEquivalentSets(Node<Key, Value>* head1, Node<Key, Value>* head2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
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