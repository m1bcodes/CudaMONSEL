#include "LinkedList.cu"

#include <stdio.h>

template<typename Key, typename Value>
void Node<Key, Value>::InsertHead(Node<Key, Value>** head, Key k, Value v)
{
   //printf("inserting %c: %f\n", k, v);
   Node<Key, Value>* newOne = new Node<Key, Value>();
   newOne->key = k;
   newOne->val = v;
   newOne->next = NULL;

   if (*head == NULL) {
      *head = newOne;
   }
   else {
      newOne->next = *head;
      *head = newOne;
   }
}

template<typename Key, typename Value>
void Node<Key, Value>::RemoveHead(Node<Key, Value>** head)
{
   if (head == NULL) {
      return;
   }

   Node<Key, Value>* tmp = (*head)->next;
   //printf("removing %c: %f\n", (*head)->key, (*head)->val);
   delete *head;
   *head = tmp;
}

template<typename Key, typename Value>
void Node<Key, Value>::RemoveAll(Node<Key, Value>** head)
{
   while (*head != NULL) {
      Node<Key, Value>::RemoveHead(head);
   }
}

template<typename Key, typename Value>
Key Node<Key, Value>::GetKey()
{
   return key;
}

template<typename Key, typename Value>
Value Node<Key, Value>::GetValue()
{
   return val;
}

template<typename Key, typename Value>
Node<Key, Value>* Node<Key, Value>::GetNext()
{
   return next;
}

//template<typename Key, typename Value>
//void Node<Key, Value>::PrintList(Node* head)
//{
//   if (head != NULL) {
//      //printf("%c: %f\n", head->key, head->val);
//      PrintList(head->next);
//   }
//}

//__host__ __device__ void BuildListTest(Node<char, float>** head)
//{
//   Node<char, float>::InsertHead(head, 'a', 0.0f);
//   Node<char, float>::InsertHead(head, 'b', 1.0f);
//   Node<char, float>::InsertHead(head, 'c', 2.0f);
//}

//__global__ void Test1()
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//
//   Node<char, float>* head = NULL;
//   printf("A\n");
//   BuildListTest(&head);
//   printf("B\n");
//   Node<char, float>::PrintList(head);
//   printf("C\n");
//   Node<char, float>::RemoveAll(&head);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//}
