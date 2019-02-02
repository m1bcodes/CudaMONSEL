#ifndef _LINKED_LIST_CUH_
#define _LINKED_LIST_CUH_

#include "cuda_runtime.h"

template<typename Key, typename Value>
class Node
{
public:
   __host__ __device__ static void InsertHead(Node** head, Key key, Value val);
   __host__ __device__ static void RemoveHead(Node** head);

   __host__ __device__ static void Remove(Node** head, Key k, bool equals(Key, Key));
   __host__ __device__ static void RemoveAll(Node** head);

   __host__ __device__ static bool AreEquivalentNodes(Node*, Node*, bool (*equalKeys)(Key, Key), bool equalValues(Value, Value));

   __host__ __device__ static void RemoveRepeatedNodes(Node**, bool(*equalKeys)(Key, Key), bool equalValues(Value, Value));
   __host__ __device__ static bool IsSet(Node*, bool(*equalKeys)(Key, Key), bool equalValues(Value, Value));
   __host__ __device__ static bool AreEquivalentSets(Node*, Node*, bool(*equalKeys)(Key, Key), bool equalValues(Value, Value));

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
__host__ __device__ void Node<Key, Value>::Remove(Node<Key, Value>** head, Key k, bool equals(Key, Key))
{
   while (*head != NULL) {
      if (equals((*head)->GetKey(), k)) {
         Node<Key, Value>::RemoveHead(head);
      }
      else {
         head = &((*head)->next);
      }
   }
}

template<typename Key, typename Value>
__host__ __device__ bool Node<Key, Value>::AreEquivalentNodes(Node<Key, Value>* head1, Node<Key, Value>* head2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   if (head1 == NULL && head2 == NULL) {
      return true;
   }

   if (head1 == NULL || head2 == NULL) {
      return false;
   }
   return (equalKeys(head1->GetKey(), head2->GetKey()) && equalValues(head1->GetValue(), head2->GetValue()));
}

template<typename Key, typename Value>
__host__ __device__ void Node<Key, Value>::RemoveRepeatedNodes(Node<Key, Value>** head, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   if (IsSet(*head, equalKeys, equalValues)) {
      return;
   }

   Node<Key, Value>** head1 = head;
   while ((*head1) != NULL) {
      Node<Key, Value>** head2 = &((*head1)->next);
      while ((*head2) == NULL) {
         if (AreEquivalentNodes(*head1, *head2, equalKeys, equalValues)) {
            RemoveHead(head2);
         }
         else {
            head2 = &((*head2)->next);
         }
      }
      head1 = &((*head1)->next);
   }
}

template<typename Key, typename Value>
__host__ __device__ bool Node<Key, Value>::IsSet(Node<Key, Value>* head, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   Node<Key, Value>* head1 = head;
   while (head1 != NULL) {
      Node<Key, Value>* head2 = head1->GetNext();
      while (head2 != NULL) {
         if (AreEquivalentNodes(head1, head2, equalKeys, equalValues)) {
            return false;
         }
         head2 = head2->GetNext();
      }
      head1 = head1->GetNext();
   }
   return true;
}

template<typename Key, typename Value>
__host__ __device__ bool Node<Key, Value>::AreEquivalentSets(Node<Key, Value>* h1, Node<Key, Value>* h2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
{
   if (!IsSet(h1, equalKeys, equalValues) || !IsSet(h2, equalKeys, equalValues)) {
      return false;
   }

   Node<Key, Value>* head1 = h1;
   while (head1 != NULL) {
      Node<Key, Value>* head2 = h2;
      while (true) {
         if (head2 == NULL) {
            return false;
         }
         if (AreEquivalentNodes(head1, head2, equalKeys, equalValues)) {
            break;
         }
         head2 = head2->GetNext();
      }
      head1 = head1->GetNext();
   }
   return true;
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

//__host__ __device__ void BuildList1(Node<String, float>** head)
//{
//   String a("a");
//   String b("b");
//   String c("c");
//   Node<String, float>::InsertHead(head, a, 0.0f);
//   Node<String, float>::InsertHead(head, b, 1.0f);
//   Node<String, float>::InsertHead(head, c, 2.0f);
//}
//
//__host__ __device__ void BuildList2(Node<String, float>** head)
//{
//   String a("a");
//   String b("b");
//   String a1("a");
//   Node<String, float>::InsertHead(head, a, 0.0f);
//   Node<String, float>::InsertHead(head, b, 1.0f);
//   Node<String, float>::InsertHead(head, a1, 0.1f);
//}
//
//__host__ __device__ void PrintListInOrder(Node<String, float>* head)
//{
//   if (head != NULL) {
//      printf("%s: %f\n", head->GetKey().Get(), head->GetValue());
//      PrintListInOrder(head->GetNext());
//   }
//}
//
//__global__ void kernel(String::pAreEqual h_strCmpPointFunction)
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//   Node<String, float>* head1 = NULL;
//   Node<String, float>* head2 = NULL;
//   printf("A\n");
//   BuildList1(&head1);
//   BuildList2(&head2);
//   printf("B\n");
//   PrintListInOrder(head1);
//   PrintListInOrder(head2);
//   printf("C\n");
//   printf("%d\n", Node<String, float>::IsSet(head1, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
//   printf("%d\n", Node<String, float>::IsSet(head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
//   printf("D\n");
//   printf("%d\n", Node<String, float>::AreEquivalentSets(head1, head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
//   Node<String, float>::RemoveRepeatedNodes(&head2, h_strCmpPointFunction, [](float a, float b) { return a == b; });
//   PrintListInOrder(head2);
//   printf("E\n");
//   Node<String, float>::Remove(&head2, String("a"), h_strCmpPointFunction);
//   PrintListInOrder(head2);
//   printf("%d\n", Node<String, float>::AreEquivalentSets(head1, head2, h_strCmpPointFunction, [](float a, float b) { return a == b; }));
//   Node<String, float>::RemoveAll(&head1);
//   Node<String, float>::RemoveAll(&head2);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//}
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

#endif