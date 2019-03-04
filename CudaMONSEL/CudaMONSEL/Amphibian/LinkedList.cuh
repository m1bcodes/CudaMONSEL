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
      //Node<T>* newOne = NULL;
      //if (*headAddr == NULL) {
      //   newOne = new Node<T>(v, NULL);
      //}
      //else {
      //   newOne = new Node<T>(v, *headAddr);
      //}
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
   __host__ __device__ bool Exists(Node<T>* head, T target, bool (*cmp)(T, T))
   {
      while (head != NULL) {
         if (cmp(head->GetValue(), target)) {
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

   template<typename T>
   __host__ __device__ void AddAllAsSet(Node<T>** res, Node<T>* other, bool(*cmp)(T, T))
   {
      while (other != NULL) {
         Node<T>* resHead = *res;
         bool found = false;
         while (resHead != NULL) {
            if (cmp(resHead->GetValue(), other->GetValue())) {
               found = true;
               break;
            }
            resHead = resHead->GetNext();
         }
         if (!found) {
            InsertHead(res, other->GetValue());
         }
         other = other->GetNext();
      }
   }

   template<typename T>
   __host__ __device__ void BuildList(Node<T>** res, T list[], int len) {
      for (int k = 0; k < len; ++k) {
         InsertHead(res, list[k]);
      }
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
         //InsertHead<KeyT, ValueT>(newHeadAddr, head->GetKey(), head->GetValue());
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
      auto r = GetValue(head, k, equalKeys);
      return !(*((int*)&r) == NULL);
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
         auto k = dataHead->GetKey();
         if (!LinkedList::Exists<K>(*headAddr, k, KeyCmp)) {
            LinkedList::InsertHead<K>(headAddr, k);
         }
         dataHead = dataHead->GetNext();
      }
   }
}


#endif