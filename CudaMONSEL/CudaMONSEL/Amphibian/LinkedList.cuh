#ifndef _LINKED_LIST_CUH_
#define _LINKED_LIST_CUH_

#include "cuda_runtime.h"

namespace LinkedList
{
   template<typename T>
   class Node
   {
   public:
      __host__ __device__ Node();
      __host__ __device__ Node(T, Node*);

      __host__ __device__ T GetValue();
      __host__ __device__ Node* GetNext();
      __host__ __device__ Node** GetNextAddr();
      __host__ __device__ void UpdateNext(Node*);

   private:
      T val;
      Node* next;
   };

   template<typename T>
   __host__ __device__ Node<T>::Node()
   {
   }

   template<typename T>
   __host__ __device__ Node<T>::Node(T v, Node* n) : val(v), next(n)
   {
   }

   template<typename T>
   __host__ __device__ T Node<T>::GetValue()
   {
      return val;
   }

   template<typename T>
   __host__ __device__ Node<T>* Node<T>::GetNext()
   {
      return next;
   }

   template<typename T>
   __host__ __device__ void Node<T>::UpdateNext(Node<T>* newNext)
   {
      next = newNext;
   }

   template<typename T>
   __host__ __device__ Node<T>** Node<T>::GetNextAddr()
   {
      return &next;
   }

   template<typename T>
   __host__ __device__ void InsertHead(Node<T>** headAddr, T v)
   {
      Node<T>* newOne = (*headAddr == NULL) ? new Node<T>(v, NULL) : new Node<T>(v, *headAddr);
      *headAddr = newOne;
   }

   template<typename T>
   __host__ __device__ void InsertNext(Node<T>** headAddr, T v)
   {
      Node<T>* newOne;
      if ((*headAddr) == NULL) {
         newOne = new Node<T>(v, NULL);
         (*headAddr) = newOne;
      }
      else {
         newOne = new Node<T>(v, (*headAddr)->GetNext());
         (*headAddr)->UpdateNext(newOne);
      }
   }

   template<typename T>
   __host__ __device__ Node<T>* DeepCopy(Node<T>* head)
   {
      Node<T>* newHead = NULL;
      Node<T>** newHeadAddr = &newHead;
      while (head != NULL) {
         InsertNext<T>(newHeadAddr, head->GetValue());
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
      return newHead;
   }

   template<typename T>
   __host__ __device__ T RemoveHead(Node<T>** headAddr)
   {
      if (*headAddr == NULL) {
         return NULL;
      }

      T v = (*headAddr)->GetValue();

      Node<T>* next = (*headAddr)->GetNext();
      delete (*headAddr);
      *headAddr = next;

      return v;
   }

   template<typename T>
   __host__ __device__ void RemoveAll(Node<T>** headAddr)
   {
      while (*headAddr != NULL) {
         RemoveHead(headAddr);
      }
   }
}

namespace LinkedListKV
{
   template<typename Key, typename Value>
   class Node
   {
   public:
      __host__ __device__ Node();
      __host__ __device__ Node(Key, Value, Node*);

      __host__ __device__ Key GetKey();
      __host__ __device__ Value GetValue();
      __host__ __device__ Node* GetNext();
      __host__ __device__ Node** GetNextAddr();
      __host__ __device__ void UpdateNext(Node* newNext);

   private:
      Key key;
      Value val;
      Node* next;
   };

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

   template<typename Key, typename Value>
   __host__ __device__ Node<Key, Value>** Node<Key, Value>::GetNextAddr()
   {
      return &next;
   }

   template<typename Key, typename Value>
   __host__ __device__ void Node<Key, Value>::UpdateNext(Node<Key, Value>* newNext)
   {
      next = newNext;
   }

   template<typename Key, typename Value>
   __host__ __device__ void InsertHead(Node<Key, Value>** head, Key k, Value v)
   {
      Node<Key, Value>* newOne = (*head == NULL) ? new Node<Key, Value>(k, v, NULL) : new Node<Key, Value>(k, v, *head);
      *head = newOne;
   }

   template<typename Key, typename Value>
   __host__ __device__ void InsertNext(Node<Key, Value>** head, Key k, Value v)
   {
      Node<Key, Value>* newOne;
      if ((*head == NULL)) {
         newOne = new Node<Key, Value>(k, v, NULL);
         (*head) = newOne;
      }
      else {
         newOne = new Node<Key, Value>(k, v, *head);
         (*head)->UpdateNext(newOne);
      }
   }

   template<typename typename Key, typename Value>
   __host__ __device__ Node<Key, Value>* DeepCopy(Node<Key, Value>* head)
   {
      Node<Key, Value>* newHead = NULL;
      Node<Key, Value>** newHeadAddr = &newHead;
      while (head != NULL) {
         InsertNext<Key, Value>(*newHeadAddr, head->GetKey(), head->GetValue());
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
      return newHead;
   }

   template<typename Key, typename Value>
   __host__ __device__ Value RemoveHead(Node<Key, Value>** head)
   {
      if (head == NULL) {
         return NULL;
      }

      Value v = (*head)->GetValue();

      Node<Key, Value>* tmp = (*head)->GetNext();
      delete *head;
      *head = tmp;

      return v;
   }

   template<typename Key, typename Value>
   __host__ __device__ Value Remove(Node<Key, Value>** head, Key k, bool equals(Key, Key))
   {
      while (*head != NULL) {
         if (equals((*head)->GetKey(), k)) {
            return RemoveHead(head);
         }
         else {
            head = (*head)->GetNextAddr();
         }
      }
      return NULL;
   }

   template<typename Key, typename Value>
   __host__ __device__ Value GetValue(Node<Key, Value>* head, Key k, bool equalKeys(Key, Key))
   {
      while (head != NULL) {
         if (equalKeys(head->GetKey(), k)) {
            return head->GetValue();
         }
         head = head->GetNext();
      }
      return NULL;
   }

   template<typename Key, typename Value>
   __host__ __device__ bool ContainsKey(Node<Key, Value>* head, Key k, bool(*equalKeys)(Key, Key))
   {
      return !(GetValue(head, k, equalKeys) == NULL);
   }

   template<typename Key, typename Value>
   __host__ __device__ bool AreEquivalentNodes(Node<Key, Value>* head1, Node<Key, Value>* head2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
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
   __host__ __device__ void RemoveRepeatedNodes(Node<Key, Value>** head, bool equalKeys(Key, Key), bool equalValues(Value, Value))
   {
      if (IsSet(*head, equalKeys, equalValues)) {
         return;
      }

      Node<Key, Value>** head1 = head;
      while ((*head1) != NULL) {
         Node<Key, Value>** head2 = (*head1)->GetNextAddr();
         while ((*head2) == NULL) {
            if (AreEquivalentNodes(*head1, *head2, equalKeys, equalValues)) {
               RemoveHead(head2);
            }
            else {
               head2 = (*head2)->GetNextAddr();
            }
         }
         head1 = (*head1)->GetNextAddr();
      }
   }

   template<typename Key, typename Value>
   __host__ __device__ bool IsSet(Node<Key, Value>* head, bool equalKeys(Key, Key), bool equalValues(Value, Value))
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
   __host__ __device__ bool AreEquivalentSets(Node<Key, Value>* h1, Node<Key, Value>* h2, bool equalKeys(Key, Key), bool equalValues(Value, Value))
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
   __host__ __device__ void RemoveAll(Node<Key, Value>** head)
   {
      while (*head != NULL) {
         RemoveHead(head);
      }
   }
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