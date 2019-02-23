#ifndef _LINKED_LIST_CUH_
#define _LINKED_LIST_CUH_

#include <cuda_runtime.h>

namespace LinkedList
{
   template<typename T>
   class Node
   {
   public:
      __host__ __device__ Node();
      __host__ __device__ Node(T, Node*);

      __host__ __device__ Node& operator=(const Node&);

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
   __host__ __device__ Node<T>& Node<T>::operator=(const Node<T>& rhs)
   {
      val = rhs.val;
      next = rhs.next;
      return *this;
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
   __host__ __device__ void DeepCopy(Node<T>** newHeadAddr, Node<T>* head)
   {
      while (head != NULL) {
         InsertNext<T>(newHeadAddr, head->GetValue());
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
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

   template<typename T>
   __host__ __device__ bool Exists(Node<T>* head, T target, bool (*KeyCmp)(T, T))
   {
      while (head != NULL) {
         if (KeyCmp(head->GetValue(), target)) {
            return true;
         }
         head = head->GetNext();
      }
      return false;
   }

   template<typename T>
   __host__ __device__ int Size(Node<T>* head)
   {
      int c = 0;
      while (head != NULL) {
         head = head->GetNext();
         c++;
      }
      return c;
   }
}

namespace LinkedListKV
{
   template<typename KeyT, typename ValueT>
   class Node
   {
   public:
      __host__ __device__ Node();
      __host__ __device__ Node(KeyT, ValueT, Node*);

      __host__ __device__ Node& operator=(const Node&);

      __host__ __device__ KeyT GetKey();
      __host__ __device__ ValueT GetValue();
      __host__ __device__ Node* GetNext();
      __host__ __device__ Node** GetNextAddr();
      __host__ __device__ void UpdateNext(Node* newNext);

   private:
      KeyT key;
      ValueT val;
      Node* next;
   };

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>::Node<KeyT, ValueT>()
   {
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>::Node<KeyT, ValueT>(KeyT k, ValueT v, Node* n) : key(k), val(v), next(n)
   {
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>& Node<KeyT, ValueT>::operator=(const Node<KeyT, ValueT>& rhs)
   {
      key = rhs.key;
      val = rhs.val;
      next = rhs.next;
   }
   
   template<typename KeyT, typename ValueT>
   __host__ __device__ KeyT Node<KeyT, ValueT>::GetKey()
   {
      return key;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ ValueT Node<KeyT, ValueT>::GetValue()
   {
      return val;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>* Node<KeyT, ValueT>::GetNext()
   {
      return next;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>** Node<KeyT, ValueT>::GetNextAddr()
   {
      return &next;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void Node<KeyT, ValueT>::UpdateNext(Node<KeyT, ValueT>* newNext)
   {
      next = newNext;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void InsertHead(Node<KeyT, ValueT>** head, KeyT k, ValueT v)
   {
      Node<KeyT, ValueT>* newOne = (*head == NULL) ? new Node<KeyT, ValueT>(k, v, NULL) : new Node<KeyT, ValueT>(k, v, *head);
      *head = newOne;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void InsertNext(Node<KeyT, ValueT>** head, KeyT k, ValueT v)
   {
      Node<KeyT, ValueT>* newOne;
      if ((*head == NULL)) {
         newOne = new Node<KeyT, ValueT>(k, v, NULL);
         (*head) = newOne;
      }
      else {
         newOne = new Node<KeyT, ValueT>(k, v, *head);
         (*head)->UpdateNext(newOne);
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void DeepCopy(Node<KeyT, ValueT>** newHeadAddr, Node<KeyT, ValueT>* head)
   {
      while (head != NULL) {
         InsertNext<KeyT, ValueT>(newHeadAddr, head->GetKey(), head->GetValue());
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ ValueT RemoveHead(Node<KeyT, ValueT>** head)
   {
      if (head == NULL) {
         return NULL;
      }

      ValueT v = (*head)->GetValue();

      Node<KeyT, ValueT>* tmp = (*head)->GetNext();
      delete *head;
      *head = tmp;

      return v;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ ValueT Remove(Node<KeyT, ValueT>** head, KeyT k, bool equals(KeyT, KeyT))
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

   template<typename KeyT, typename ValueT>
   __host__ __device__ ValueT GetValue(Node<KeyT, ValueT>* head, KeyT target, bool equalKeys(KeyT, KeyT))
   {
      while (head != NULL) {
         if (equalKeys(head->GetKey(), target)) {
            return head->GetValue();
         }
         head = head->GetNext();
      }
      return NULL;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool ContainsKey(Node<KeyT, ValueT>* head, KeyT k, bool(*equalKeys)(KeyT, KeyT))
   {
      return !(GetValue(head, k, equalKeys) == NULL);
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool AreEquivalentNodes(Node<KeyT, ValueT>* head1, Node<KeyT, ValueT>* head2, bool equalKeys(KeyT, KeyT), bool equalValues(ValueT, ValueT))
   {
      if (head1 == NULL && head2 == NULL) {
         return true;
      }

      if (head1 == NULL || head2 == NULL) {
         return false;
      }
      return (equalKeys(head1->GetKey(), head2->GetKey()) && equalValues(head1->GetValue(), head2->GetValue()));
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void RemoveRepeatedNodes(Node<KeyT, ValueT>** head, bool equalKeys(KeyT, KeyT), bool equalValues(ValueT, ValueT))
   {
      if (IsSet(*head, equalKeys, equalValues)) {
         return;
      }

      Node<KeyT, ValueT>** head1 = head;
      while ((*head1) != NULL) {
         Node<KeyT, ValueT>** head2 = (*head1)->GetNextAddr();
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

   template<typename KeyT, typename ValueT>
   __host__ __device__ void RemoveAll(Node<KeyT, ValueT>** head)
   {
      while (*head != NULL) {
         RemoveHead(head);
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ int Size(Node<KeyT, ValueT>* head)
   {
      int sz = 0;
      while (head != NULL) {
         head = head->GetNext();
         sz++;
      }
      return sz;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool IsSet(Node<KeyT, ValueT>* head, bool equalKeys(KeyT, KeyT), bool equalValues(ValueT, ValueT))
   {
      Node<KeyT, ValueT>* head1 = head;
      while (head1 != NULL) {
         Node<KeyT, ValueT>* head2 = head1->GetNext();
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

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool AreEquivalentSets(Node<KeyT, ValueT>* h1, Node<KeyT, ValueT>* h2, bool equalKeys(KeyT, KeyT), bool equalValues(ValueT, ValueT))
   {
      if (!IsSet(h1, equalKeys, equalValues) || !IsSet(h2, equalKeys, equalValues)) {
         return false;
      }

      Node<KeyT, ValueT>* head1 = h1;
      while (head1 != NULL) {
         Node<KeyT, ValueT>* head2 = h2;
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
}

// Advanced Templates

namespace AdvancedLinkedList
{
   template<typename K, typename V>
   __host__ __device__ void AddAllKeys(LinkedList::Node<K>** headAddr, LinkedListKV::Node<K, V>* dataHead, bool (*KeyCmp)(K, K))
   {
      if (dataHead == NULL) {
         return;
      }
      while (dataHead != NULL) {
         if (!LinkedList::Exists<K>(*headAddr, dataHead->GetKey(), KeyCmp)) {
            LinkedList::InsertHead<K>(headAddr, dataHead->GetKey());
         }
         dataHead = dataHead->GetNext();
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

//__device__ void PrintList(LinkedList::Node<int>* head)
//{
//   while (head != NULL) {
//      printf("- %d\n", head->GetValue());
//      head = head->GetNext();
//   }
//   printf("-------------------\n");
//}
//
//__device__ void doThings()
//{
//   LinkedList::Node<int> * head = NULL;
//   LinkedList::Node<int>** newHeadAddr = &head;
//   LinkedList::InsertNext(newHeadAddr, 0);
//   LinkedList::InsertNext(newHeadAddr, 1);
//   newHeadAddr = (*newHeadAddr)->GetNextAddr();
//   LinkedList::InsertNext(newHeadAddr, 2);
//
//   PrintList(head);
//   LinkedList::Node<int> * head2 = LinkedList::DeepCopy(head);
//   LinkedList::RemoveAll(&head);
//   PrintList(head);
//   PrintList(head2);
//   LinkedList::RemoveAll(&head2);
//   PrintList(head2);
//}

#endif