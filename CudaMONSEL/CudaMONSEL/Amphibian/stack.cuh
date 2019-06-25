#ifndef _STACK_CUH_
#define _STACK_CUH_

#include <cuda_runtime.h>

#include "Amphibian\LinkedList.cuh"

namespace amp
{
   template<typename T>
   class stack
   {
   public:
      __host__ __device__ stack();
      __host__ __device__ stack(const stack&);
      __host__ __device__ ~stack();

      __host__ __device__ const T& top() const;

      __host__ __device__ void pop();
      __host__ __device__ void push(const T& t);

      //__host__ __device__ unsigned int size() const;
      __host__ __device__ bool empty() const;

   private:
      __host__ __device__ void DeepCopy(const stack&);
      __host__ __device__ void RemoveAll();

      DLinkedList::Node<T>* head;
   };

   template<typename T>
   __host__ __device__ stack<T>::stack() : head(nullptr)
   {
   }

   template<typename T>
   __host__ __device__ stack<T>::stack(const stack& other) : head(nullptr)
   {
      DeepCopy(other);
   }

   template<typename T>
   __host__ __device__ stack<T>::~stack()
   {
      RemoveAll();
   }

   template<typename T>
   __host__ __device__ const T& stack<T>::top() const
   {
      return head->GetValue();
   }

   template<typename T>
   __host__ __device__ void stack<T>::pop()
   {
      DLinkedList::Node<T>::Remove(&head);
   }

   template<typename T>
   __host__ __device__ void stack<T>::push(const T& t)
   {
      DLinkedList::Node<T>::Insert(&head, t);
   }

   //template<typename T>
   //__host__ __device__ unsigned int stack<T>::size()
   //{
   //   return sz;
   //}

   template<typename T>
   __host__ __device__ bool stack<T>::empty() const
   {
      return head == nullptr;
   }

   template<typename T>
   __host__ __device__ void stack<T>::DeepCopy(const stack& other)
   {
      DLinkedList::Node<T>* itr = other.head;
      while (itr) {
         DLinkedList::Node<T>::Insert(&head, itr->val);
         itr = itr->GetNext();
      }
   }

   template<typename T>
   __host__ __device__ void stack<T>::RemoveAll()
   {
      while (head) {
         DLinkedList::Node<T>::Remove(&head);
      }
   }
}

#endif