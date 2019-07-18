#ifndef _LINKED_LIST_CUH_
#define _LINKED_LIST_CUH_

#include "Amphibian/Hasher.cuh"

#include <cuda_runtime.h>

#include <stdio.h>

namespace LinkedList
{
   template<typename T>
   class Node
   {
   public:
      //__host__ __device__ Node();
      __host__ __device__ Node(const T&, Node*);

      __host__ __device__ Node& operator=(const Node&);

      __host__ __device__ T& GetValue();
      __host__ __device__ const T& GetValue() const;
      __host__ __device__ Node* GetNext();
      __host__ __device__ const Node* GetNext() const;
      __host__ __device__ Node** GetNextAddr();
      __host__ __device__ void UpdateNext(Node*);

   private:
      T val;
      Node* next;
   };

   //template<typename T>
   //__host__ __device__ Node<T>::Node()
   //{
   //}

   template<typename T>
   __host__ __device__ Node<T>::Node(const T& v, Node* n) : val(v), next(n)
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
   __host__ __device__ T& Node<T>::GetValue()
   {
      return val;
   }

   template<typename T>
   __host__ __device__ const T& Node<T>::GetValue() const
   {
      return val;
   }

   template<typename T>
   __host__ __device__ Node<T>* Node<T>::GetNext()
   {
      return next;
   }

   template<typename T>
   __host__ __device__ const Node<T>* Node<T>::GetNext() const
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
   __host__ __device__ void InsertHead(Node<T>** headAddr, const T& v)
   {
      Node<T>* newOne = (*headAddr == nullptr) ? new Node<T>(v, nullptr) : new Node<T>(v, *headAddr);
      if (!newOne) printf("LinkedList::InsertHead: overflowed\n");
      *headAddr = newOne;
   }

   template<typename T>
   __host__ __device__ void InsertNext(Node<T>** headAddr, T& v)
   {
      Node<T>* newOne;
      if ((*headAddr) == nullptr) {
         newOne = new Node<T>(v, nullptr);
         if (!newOne) printf("LinkedList::InsertNext 1: overflowed\n");
         (*headAddr) = newOne;
      }
      else {
         newOne = new Node<T>(v, (*headAddr)->GetNext());
         if (!newOne) printf("LinkedList::InsertNext 2: overflowed\n");
         (*headAddr)->UpdateNext(newOne);
      }
   }

   template<typename T>
   __host__ __device__ void DeepCopy(Node<T>** newHeadAddr, Node<T>* head)
   {
      while (head != nullptr) {
         T v = head->GetValue();
         InsertNext<T>(newHeadAddr, v);
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
   }

   template<typename T>
   __host__ __device__ bool RemoveHead(Node<T>** headAddr, T& ret)
   {
      if (*headAddr == nullptr) {
         return false;
      }

      ret = (*headAddr)->GetValue();

      Node<T>* next = (*headAddr)->GetNext();
      delete (*headAddr);
      *headAddr = next;

      return true;
   }

   template<typename T>
   __host__ __device__ bool RemoveHead(Node<T>** headAddr)
   {
      if (*headAddr == nullptr) {
         return false;
      }

      Node<T>* next = (*headAddr)->GetNext();
      delete (*headAddr);
      *headAddr = next;

      return true;
   }

   template<typename T>
   __host__ __device__ void RemoveAll(Node<T>** headAddr)
   {
      while (*headAddr != nullptr) {
         RemoveHead(headAddr);
      }
   }

   template<typename T, typename TCompare>
   __host__ __device__ T Remove(Node<T>** head, const T& k)
   {
      TCompare cmp;
      while (*head != nullptr) {
         T tmpK = (*head)->GetValue();
         if (cmp(tmpK, k)) {
            T ret;
            RemoveHead(head, ret);
            return ret;
         }
         else {
            head = (*head)->GetNextAddr();
         }
      }
      return NULL;
   }

   template<typename T, typename TCompare>
   __host__ __device__ bool Exists(Node<T> * head, const T& target)
   {
      TCompare cmp;
      while (head != nullptr) {
         T tmpV = head->GetValue();
         if (cmp(tmpV, target)) {
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
      while (head != nullptr) {
         head = head->GetNext();
         c++;
      }
      return c;
   }

   template<typename T>
   __host__ __device__ void AddAllAsSet(Node<T>** res, Node<T>* other, bool(*cmp)(T&, T&))
   {
      while (other != nullptr) {
         Node<T>* resHead = *res;
         bool found = false;
         while (resHead != nullptr) {
            T tmpV1 = resHead->GetValue();
            T tmpV2 = other->GetValue();
            if (cmp(tmpV1, tmpV2)) {
               found = true;
               break;
            }
            resHead = resHead->GetNext();
         }
         if (!found) {
            T v = other->GetValue();
            InsertHead(res, v);
         }
         other = other->GetNext();
      }
   }

   template<typename T>
   __host__ __device__ void BuildList(Node<T>** res, T list[], int len)
   {
      for (int k = 0; k < len; ++k) {
         InsertHead(res, list[k]);
      }
   }

   template<typename T, typename THasher>
   __host__ __device__ unsigned int HashCode(Node<T>* list)
   {
      THasher hasher;
      unsigned int res = 0;
      while (list != nullptr) {
         T v = list->GetValue();
         res += hasher(v);
         list = list->GetNext();
      }
      return res;
   }
}

namespace LinkedListKV
{
   template<typename KeyT, typename ValueT>
   class Node
   {
   public:
      __host__ __device__ Node();
      __host__ __device__ Node(const KeyT&, const ValueT&, Node*);

      __host__ __device__ Node& operator=(const Node&);

      __host__ __device__ KeyT& GetKey();
      __host__ __device__ ValueT& GetValue();
      __host__ __device__ const KeyT& GetKey() const;
      __host__ __device__ const ValueT& GetValue() const;
      __host__ __device__ Node* GetNext();
      __host__ __device__ const Node* GetNext() const;
      __host__ __device__ Node** GetNextAddr();
      __host__ __device__ void UpdateNext(Node* newNext);
      __host__ __device__ void MapVal(const ValueT&, const ValueT&(*mapper)(const ValueT&, ValueT&));

      class First
      {
      public:
         __host__ __device__ First(Node& node) : node(node) {}
         __host__ __device__ operator KeyT&() const { return (KeyT&)node.key; }

      private:
         Node& node;
      } first;

      class Second
      {
      public:
         __host__ __device__ Second(Node& node) : node(node) {}
         __host__ __device__ operator ValueT&() const { return (ValueT&)node.val; }

      private:
         Node& node;
      } second;

   private:
      KeyT key;
      ValueT val;
      Node* next;
   };

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>::Node<KeyT, ValueT>() : first(*this), second(*this)
   {
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>::Node<KeyT, ValueT>(const KeyT& k, const ValueT& v, Node* n) : key(k), val(v), next(n), first(*this), second(*this)
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
   __host__ __device__ KeyT& Node<KeyT, ValueT>::GetKey()
   {
      return key;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ ValueT& Node<KeyT, ValueT>::GetValue()
   {
      return val;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ const KeyT& Node<KeyT, ValueT>::GetKey() const
   {
      return key;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ const ValueT& Node<KeyT, ValueT>::GetValue() const
   {
      return val;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ Node<KeyT, ValueT>* Node<KeyT, ValueT>::GetNext()
   {
      return next;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ const Node<KeyT, ValueT>* Node<KeyT, ValueT>::GetNext() const
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

   // ValueT needs to have copy constructor
   template<typename KeyT, typename ValueT>
   __host__ __device__ void Node<KeyT, ValueT>::MapVal(const ValueT& v, const ValueT&(*mapper)(const ValueT&, ValueT&))
   {
      val = mapper(v, val);
   }

   template<typename KeyT, typename ValueT, typename KHasher, typename VHasher>
   __host__ __device__ unsigned int HashCode(Node<KeyT, ValueT>* list)
   {
      KHasher khasher;
      VHasher vhasher;
      unsigned int res = 0;
      while (list != nullptr) {
         KeyT k = list->GetKey();
         ValueT v = list->GetValue();
         unsigned int h = khasher(k) ^ vhasher(v); // https://docs.oracle.com/javase/7/docs/api/java/util/Map.Entry.html#hashCode()
         res += h; // https://docs.oracle.com/javase/7/docs/api/java/util/Map.html#hashCode()
         list = list->GetNext();
      }
      return res;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void InsertHead(Node<KeyT, ValueT>** head, const KeyT& k, const ValueT& v)
   {
      Node<KeyT, ValueT>* newOne = (*head == nullptr) ? new Node<KeyT, ValueT>(k, v, nullptr) : new Node<KeyT, ValueT>(k, v, *head);
      if (!newOne) printf("LinkedListKV::InsertHead: overflowed\n");
      *head = newOne;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void InsertHead(Node<KeyT, ValueT>** head, Node<KeyT, ValueT>* newOne)
   {
      newOne->UpdateNext(*head);
      *head = newOne;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void InsertNext(Node<KeyT, ValueT>** head, KeyT& k, ValueT& v)
   {
      Node<KeyT, ValueT>* newOne;
      if ((*head == nullptr)) {
         newOne = new Node<KeyT, ValueT>(k, v, nullptr);
         if (!newOne) printf("LinkedListKV::InsertNext 1: overflowed\n");
         (*head) = newOne;
      }
      else {
         newOne = new Node<KeyT, ValueT>(k, v, *head);
         if (!newOne) printf("LinkedListKV::InsertNext 2: overflowed\n");
         (*head)->UpdateNext(newOne);
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void DeepCopy(Node<KeyT, ValueT>** newHeadAddr, Node<KeyT, ValueT>* head)
   {
      while (head != nullptr) {
         InsertNext<KeyT, ValueT>(newHeadAddr, head->GetKey(), head->GetValue());
         newHeadAddr = (*newHeadAddr)->GetNextAddr();
         head = head->GetNext();
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool RemoveHead(Node<KeyT, ValueT>** head, ValueT& ret)
   {
      if (*head == nullptr) {
         return false;
      }

      ret = (*head)->GetValue();

      Node<KeyT, ValueT>* tmp = (*head)->GetNext();
      delete *head;
      *head = tmp;

      return true;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ bool RemoveHead(Node<KeyT, ValueT>** head)
   {
      ValueT v;
      return RemoveHead(head, v);
   }

   template<typename KeyT, typename ValueT, typename KCompare>
   __host__ __device__ bool Remove(Node<KeyT, ValueT>** head, const KeyT& k, ValueT& ret)
   {
      KCompare kcmp;
      while (*head != nullptr) {
         KeyT tmpK = (*head)->GetKey();
         if (kcmp(tmpK, k)) {
            ValueT ret;
            RemoveHead(head, ret);
            return true;
         }
         else {
            head = (*head)->GetNextAddr();
         }
      }
      return false;
   }

   template<typename KeyT, typename ValueT, typename KCompare>
   __host__ __device__ bool GetValue(Node<KeyT, ValueT>* head, const KeyT& target, ValueT& ret)
   {
      KCompare kcmp;
      while (head != nullptr) {
         KeyT k = head->GetKey();
         if (kcmp(k, target)) {
            ret = head->GetValue();
            return true;
         }
         head = head->GetNext();
      }
      return false;
   }

   template<typename KeyT, typename ValueT, typename KCompare>
   __host__ __device__ bool ContainsKey(Node<KeyT, ValueT>* head, const KeyT& k)
   {
      ValueT v;
      return GetValue<KeyT, ValueT, KCompare>(head, k, v);
   }

   template<typename KeyT, typename ValueT, typename KCompare>
   __host__ __device__ Node<KeyT, ValueT>* Find(Node<KeyT, ValueT>* head, const KeyT& target)
   {
      KCompare kcmp;
      while (head != nullptr) {
         KeyT k = head->GetKey();
         if (kcmp(k, target)) {
            return head;
         }
         head = head->GetNext();
      }
      return nullptr;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void RemoveAll(Node<KeyT, ValueT>** head)
   {
      while (*head != nullptr) {
         RemoveHead(head);
      }
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ int Size(Node<KeyT, ValueT>* head)
   {
      int sz = 0;
      while (head != nullptr) {
         head = head->GetNext();
         ++sz;
      }
      return sz;
   }

   template<typename KeyT, typename ValueT>
   __host__ __device__ void MapVal(const ValueT& v, Node<KeyT, ValueT>* head, void(*mapper)(const ValueT&, ValueT&))
   {
      while (head != nullptr) {
         head->MapVal(v, mapper);
         head = head->GetNext();
      }
   }
}

namespace DLinkedList
{
   template<typename T>
   class Node
   {
   public:
      __host__ __device__ Node(const Node&);
      __host__ __device__ Node(Node* fwd, Node* bwd, const T&);

      __host__ __device__ const T& GetValue() const;
      __host__ __device__ Node* GetNext() const;
      __host__ __device__ Node* GetPrev() const;

      __host__ __device__ static void Insert(Node** node, Node* newOne);
      __host__ __device__ static void Insert(Node** node, const T& val);
      __host__ __device__ static void Remove(Node** node);

   private:
      Node *fwd, *bwd;
      T val;
   };

   template<typename T>
   __host__ __device__ Node<T>::Node(const Node& other) : fwd(other.fwd), bwd(other.bwd), val(other.val)
   {
   }

   template<typename T>
   __host__ __device__ Node<T>::Node(Node* fwd, Node* bwd, const T& val) : fwd(fwd), bwd(bwd), val(val)
   {
   }


   template<typename T>
   __host__ __device__ const T& Node<T>::GetValue() const
   {
      return val;
   }

   template<typename T>
   __host__ __device__ Node<T>* Node<T>::GetNext() const
   {
      return fwd;
   }

   template<typename T>
   __host__ __device__ Node<T>* Node<T>::GetPrev() const
   {
      return bwd;
   }

   template<typename T>
   __host__ __device__ void Node<T>::Insert(Node** node, const T& val)
   {
      Node *newOne = new Node(nullptr, nullptr, val);
      if (!newOne) printf("DLinkedList::Insert: overflowed\n");
      if (!(*node)) {
         *node = newOne;
         return;
      }

      newOne->fwd = *node;
      newOne->bwd = (*node)->bwd;
      (*node)->bwd = newOne;

      if (newOne->bwd) {
         (newOne->bwd)->fwd = newOne;
      }

      *node = newOne;
   }

   template<typename T>
   __host__ __device__ void Node<T>::Insert(Node** node, Node* newOne)
   {
      if (!newOne) return;
      if (!(*node)) {
         *node = newOne;
         return;
      }

      newOne->fwd = *node;
      newOne->bwd = (*node)->bwd;
      (*node)->bwd = newOne;

      if (newOne->bwd) {
         (newOne->bwd)->fwd = newOne;
      }

      *node = newOne;
   }

   template<typename T>
   __host__ __device__ void Node<T>::Remove(Node** node)
   {
      if (!(*node)) return;
      Node *f = (*node)->fwd, *b = (*node)->bwd;

      if ((*node)->fwd) ((*node)->fwd)->bwd = b;
      if ((*node)->bwd) ((*node)->bwd)->fwd = f;
      delete *node;

      *node = f ? f : b;
   }
}

namespace amp
{
   template<typename First, typename Second>
   __host__ __device__ LinkedListKV::Node<First, Second>* make_pair(const First& first, const Second& second)
   {
      LinkedListKV::Node<First, Second>* ptr = new LinkedListKV::Node<First, Second>(first, second, nullptr);
      if (!ptr) printf("LinkedListKV::Node<First, Second>* make_pair: overflowed\n");
      return ptr;
   }
}

#endif
