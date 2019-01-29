#include <stdio.h>

#include "cuda_runtime.h"

#include "gov/nist/microanalysis/NISTMonte/MonteCarloSS.cu"
//#include "gov/nist/microanalysis/NISTMonte/Electron.cu"
//#include "gov/nist/microanalysis/Utility/CSVReader.h"

#include "Amphibian\String.cuh"
#include "Amphibian\LinkedList.cuh"

#include "CudaUtil.h"
#include "ImageUtil.h"
#include "TimeUtil.h"

__global__ void spawnElectron(unsigned int *d_arr, int idx_x, int idx_y, size_t size_x, size_t size_y)
{
   int idx = idx_y * size_x + idx_x;
   //MonteCarloSS::RegionBase e(idx);
   //d_arr[idx] = e.GetId();

   d_arr[idx] = idx;
}

__global__ void spawnElectrons(unsigned int *d_arr, size_t size_x, size_t size_y)
{
   int idx_x = threadIdx.x + blockDim.x * blockIdx.x;
   int idx_y = threadIdx.y + blockDim.y * blockIdx.y;

   printf("%d, %d", idx_x, idx_y);

   spawnElectron << <1, 1 >> >(d_arr, idx_x, idx_y, size_x, size_y);
}

void PrintArray2D(unsigned int *h_arr, size_t img_x, size_t img_y)
{
   for (int k = 0; k < img_y; ++k) {
      for (int l = 0; l < img_x; ++l) {
         std::cout << h_arr[k*img_x + l] << " ";
      }
      std::cout << std::endl;
   }
   std::cout << std::endl;
}

//int main()
//{
//   checkCudaErrors(cudaSetDevice(0));
//
//   const size_t img_x = 256;
//   const size_t img_y = 256;
//   const dim3 blockSize(2, 2, 1);
//   const size_t grid_x = (img_x / blockSize.x) + ((img_x % blockSize.x > 0) ? 1 : 0);
//   const size_t grid_y = (img_y / blockSize.y) + ((img_y % blockSize.y > 0) ? 1 : 0);
//   const dim3 gridSize(grid_x, grid_y, 1);
//
//   unsigned int *d_arr;
//   checkCudaErrors(cudaMalloc((void **)&d_arr, sizeof(int) * img_x * img_y));
//
//   spawnElectrons << <gridSize, blockSize >> >(d_arr, img_x, img_y);
//   checkCudaErrors(cudaDeviceSynchronize());
//   checkCudaErrors(cudaGetLastError());
//
//   unsigned int *h_arr = new unsigned int[img_x * img_y];
//   checkCudaErrors(cudaMemcpy(h_arr, d_arr, sizeof(int) * img_x * img_y, cudaMemcpyDeviceToHost));
//   checkCudaErrors(cudaFree(d_arr));
//
//   saveResults("a.bmp", h_arr, img_x, img_y);
//   delete[] h_arr;
//
//   return 0;
//}

//__host__ __device__ void BuildListTest(Node<String, float>** head)
//{
//   String a("a");
//   String b("b");
//   String c("c");
//   Node<String, float>::InsertHead(head, a, 0.0f);
//   Node<String, float>::InsertHead(head, b, 1.0f);
//   Node<String, float>::InsertHead(head, c, 2.0f);
//}
//
//__host__ __device__ void PrintList(Node<String, float>* head)
//{
//   if (head != NULL) {
//      printf("%s: %f\n", head->GetKey().Get(), head->GetValue());
//      //PrintList(head->GetNext());
//   }
//}
//
//__global__ void Test1()
//{
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//
//   Node<String, float>* head = NULL;
//   printf("A\n");
//   BuildListTest(&head);
//   printf("B\n");
//   PrintList(head);
//   printf("C\n");
//   Node<String, float>::RemoveAll(&head);
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//   __syncthreads();
//#endif
//}

//template<typename Key, typename Value>
//class Node
//{
//public:
//   __host__ __device__ static void InsertHead(Node** head, Key key, Value val);
//   __host__ __device__ static void RemoveHead(Node** head);
//
//   __host__ __device__ static void RemoveAll(Node** head);
//
//   __host__ __device__ Node(Key, Value, Node *);
//
//   __host__ __device__ Key GetKey();
//   __host__ __device__ Value GetValue();
//   __host__ __device__ Node<Key, Value>* GetNext();
//
//private:
//   Key key;
//   Value val;
//   Node* next;
//};
//
//template<typename Key, typename Value>
//__host__ __device__ void Node<Key, Value>::InsertHead(Node<Key, Value>** head, Key k, Value v)
//{
//   //printf("inserting %c: %f\n", k, v);
//   Node<Key, Value>* newOne = new Node<Key, Value>(k, v, NULL);
//   printf("%s?\n", newOne->GetKey().Get());
//
//   if (*head == NULL) {
//      *head = newOne;
//   }
//   else {
//      newOne->next = *head;
//      *head = newOne;
//   }
//   printf("%s\n", (*head)->GetKey().Get());
//}
//
//template<typename Key, typename Value>
//__host__ __device__ void Node<Key, Value>::RemoveHead(Node<Key, Value>** head)
//{
//   if (head == NULL) {
//      return;
//   }
//
//   Node<Key, Value>* tmp = (*head)->next;
//   //printf("removing %c: %f\n", (*head)->key, (*head)->val);
//   delete *head;
//   *head = tmp;
//}
//
//template<typename Key, typename Value>
//__host__ __device__ void Node<Key, Value>::RemoveAll(Node<Key, Value>** head)
//{
//   while (*head != NULL) {
//      Node<Key, Value>::RemoveHead(head);
//   }
//}
//
//template<typename Key, typename Value>
//__host__ __device__ Node<Key, Value>::Node<Key, Value>(Key k, Value v, Node* n)
//{
//   key = k;
//   val = v;
//   next = n;
//   printf("%s\n", key.Get());
//}
//
//template<typename Key, typename Value>
//__host__ __device__ Key Node<Key, Value>::GetKey()
//{
//   return key;
//}
//
//template<typename Key, typename Value>
//__host__ __device__ Value Node<Key, Value>::GetValue()
//{
//   return val;
//}
//
//template<typename Key, typename Value>
//__host__ __device__ Node<Key, Value>* Node<Key, Value>::GetNext()
//{
//   return next;
//}

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
__host__ __device__ Node<Key, Value>::Node<Key, Value>(Key k, Value v, Node* n)
{
   key = k;
   val = v;
   next = n;
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

__global__ void Test1()
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __syncthreads();
#endif
   Node<String, float>* head = NULL;
   String a("a");
   Node<String, float>::InsertHead(&head, a, 0.0f);
   printf("%s\n", head->GetKey().Get());
   printf("%s: %f\n", head->GetKey().Get(), head->GetValue());
   Node<String, float>::RemoveAll(&head);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __syncthreads();
#endif
}

//__global__ void kernel()
//{
//   String s("aa\n");
//   printf("%s", s.Get());
//}

int main()
{
   //kernel << <1, 1 >> >();
   Test1 << < 1, 1 >> > ();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   return 0;
}
