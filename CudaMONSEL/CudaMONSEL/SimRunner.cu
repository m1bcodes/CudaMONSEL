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

__host__ __device__ void BuildListTest(Node<String, float>** head)
{
   String a("a");
   String b("b");
   String c("c");
   Node<String, float>::InsertHead(head, a, 0.0f);
   Node<String, float>::InsertHead(head, b, 1.0f);
   Node<String, float>::InsertHead(head, c, 2.0f);
}

__host__ __device__ void PrintListInOrder(Node<String, float>* head)
{
   if (head != NULL) {
      printf("%s: %f\n", head->GetKey().Get(), head->GetValue());
      PrintListInOrder(head->GetNext());
   }
}

__global__ void Test1()
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __syncthreads();
#endif

   Node<String, float>* head = NULL;
   printf("A\n");
   BuildListTest(&head);
   printf("B\n");
   PrintListInOrder(head);
   printf("C\n");
   Node<String, float>::RemoveAll(&head);

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __syncthreads();
#endif
}

__global__ void kernel(int n)
{
   char n_a[16] = "\0";
   String::IToA(n_a, n);
   printf("%s\n", n_a);
}

int main()
{
   kernel << <1, 1 >> >(12345678);
   //Test1 << < 1, 1 >> > ();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   return 0;
}
