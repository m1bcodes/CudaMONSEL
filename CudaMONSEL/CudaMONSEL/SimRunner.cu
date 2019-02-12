#include <stdio.h>

#include <cuda_runtime.h>

#include "gov/nist/microanalysis/NISTMonte/MonteCarloSS.cu"
//#include "gov/nist/microanalysis/NISTMonte/Electron.cu"
//#include "gov/nist/microanalysis/Utility/CSVReader.h"
#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"

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

__global__ void kernel()
{
   UncertainValue2::UncertainValue2 v0(0, "abc", 5);
   UncertainValue2::UncertainValue2 v1(1);
   UncertainValue2::UncertainValue2 v2(2, 10);
   UncertainValue2::UncertainValue2 v3(2, 10);

   printf("%d\n", v1.equals(&v2));
}

int main()
{
   kernel << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   return 0;
}
