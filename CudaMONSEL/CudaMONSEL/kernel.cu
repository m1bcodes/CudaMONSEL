//#include "cuda_runtime.h"
//#include "utils.h"
//
//#include <stdio.h>
//
//class MonteCarloSS
//{
//protected:
//    int i;
//public:
//    class RegionBase;
//    __host__ __device__ MonteCarloSS() { i = 0; }
//    __host__ __device__ ~MonteCarloSS();
//};
//
//class MonteCarloSS::RegionBase
//{
//protected:
//    int k;
//public:
//    __host__ __device__ RegionBase() { k = 1; }
//    __host__ __device__ ~RegionBase();
//};
//
//__device__ int add(int x, int y)
//{
//    return x + y;
//}
//
//__global__ void addKernel(int *c, const int *a, const int *b, size_t size_x, size_t size_y)
//{
//    int idx_x = threadIdx.x + blockDim.x * blockIdx.x;
//    int idx_y = threadIdx.y + blockDim.y * blockIdx.y;
//    int idx = idx_y * size_x + idx_x;
//    c[idx] = add(a[idx], b[idx]);
//}
//
//__global__ void doThingWithCuda2(unsigned int *d_arr, size_t size_x, size_t size_y)
//{
//   int idx_x = threadIdx.x + blockDim.x * blockIdx.x;
//   int idx_y = threadIdx.y + blockDim.y * blockIdx.y;
//   int idx = idx_y * size_x + idx_x;
//   //d_arr[idx] = (idx_x + 100 * idx_y);
//   //d_arr[idx] += 1;
//   atomicAdd(d_arr + idx, 1);
//}
//
//__global__ void doThingWithCuda(unsigned int *d_arr, size_t size_x, size_t size_y)
//{
//   const size_t img_x = 32;
//   const size_t img_y = 32;
//   //const size_t maxThreadsPerBlock = 1024;
//   const dim3 blockSize(16, 16, 1);
//   const size_t grid_x = (img_x / blockSize.x) + ((img_x % blockSize.x > 0) ? 1 : 0);
//   const size_t grid_y = (img_y / blockSize.y) + ((img_y % blockSize.y > 0) ? 1 : 0);
//   const dim3 gridSize(grid_x, grid_y, 1);
//
//   doThingWithCuda2<<<gridSize, blockSize>>>(d_arr, size_x, size_y);
//    //int idx_x = threadIdx.x + blockDim.x * blockIdx.x;
//    //int idx_y = threadIdx.y + blockDim.y * blockIdx.y;
//    //int idx = idx_y * size_x + idx_x;
//    //d_arr[idx] += idx_x + 100 * idx_y;
//}
//
//__global__ void squareWithCuda2d(unsigned int **d_arr, int arrLength)
//{
//    int idInArrayx = threadIdx.x + blockDim.x * blockIdx.x; // thread index
//    int idInArrayy = threadIdx.y + blockDim.y * blockIdx.y; // thread index
//
//    d_arr[idInArrayx][idInArrayy] = idInArrayx * idInArrayy;
//}
//
//// Helper function for using CUDA to add vectors in parallel.
//cudaError_t addWithCuda(int *c, const int *a, const int *b, size_t size_x, size_t size_y)
//{
//    int *dev_a = 0;
//    int *dev_b = 0;
//    int *dev_c = 0;
//    int length = size_x * size_y;
//    cudaError_t cudaStatus;
//
//    // Choose which GPU to run on, change this on a multi-GPU system.
//    cudaStatus = cudaSetDevice(0);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
//        goto Error;
//    }
//
//    // Allocate GPU buffers for three vectors (two input, one output)    .
//    cudaStatus = cudaMalloc((void**)&dev_c, length * sizeof(int));
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMalloc failed!");
//        goto Error;
//    }
//
//    cudaStatus = cudaMalloc((void**)&dev_a, length * sizeof(int));
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMalloc failed!");
//        goto Error;
//    }
//
//    cudaStatus = cudaMalloc((void**)&dev_b, length * sizeof(int));
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMalloc failed!");
//        goto Error;
//    }
//
//    // Copy input vectors from host memory to GPU buffers.
//    cudaStatus = cudaMemcpy(dev_a, a, length * sizeof(int), cudaMemcpyHostToDevice);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMemcpy failed!");
//        goto Error;
//    }
//
//    cudaStatus = cudaMemcpy(dev_b, b, length * sizeof(int), cudaMemcpyHostToDevice);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMemcpy failed!");
//        goto Error;
//    }
//
//    // Launch a kernel on the GPU with one thread for each element.
//    const size_t maxThreadsPerBlock = 1024;
//    const dim3 blockSize(2, 2, 1);
//    const dim3 gridSize((size_x + blockSize.x - 1) / blockSize.x, (size_y + blockSize.y - 1) / blockSize.y, 1);
//    addKernel << <gridSize, blockSize >> >(dev_c, dev_a, dev_b, size_x, size_y);
//
//    // Check for any errors launching the kernel
//    cudaStatus = cudaGetLastError();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//        goto Error;
//    }
//
//    // cudaDeviceSynchronize waits for the kernel to finish, and returns
//    // any errors encountered during the launch.
//    cudaStatus = cudaDeviceSynchronize();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
//        goto Error;
//    }
//
//    // Copy output vector from GPU buffer to host memory.
//    cudaStatus = cudaMemcpy(c, dev_c, length * sizeof(int), cudaMemcpyDeviceToHost);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaMemcpy failed!");
//        goto Error;
//    }
//
//Error:
//    cudaFree(dev_c);
//    cudaFree(dev_a);
//    cudaFree(dev_b);
//
//    return cudaStatus;
//}
//
////__host__ __device__
////float gaussian1d(float x, float sigma)
////{
////    float PI = 3.141526;
////    float variance = pow(sigma, 2);
////    float exponent = -pow(x, 2) / (2 * variance);
////    return expf(exponent) / sqrt(2 * PI * variance);
////}
////
////inline __device__
////float gaussian1d_gpu(float x, float sigma)
////{
////    float PI = 3.141526;
////    float variance = __powf(sigma, 2);
////    //this doesnt work for some reason: __powf(x,2.0f)
////    float power = pow(x, 2);
////    float exponent = -power / (2 * variance);
////    return __expf(exponent) / sqrt(2 * PI * variance);
////}
////__host__ __device__
////float gaussian2d(float x, float y, float sigma)
////{
////    float PI = 3.141526;
////    float variance = pow(sigma, 2);
////    float exponent = -(pow(x, 2) + pow(y, 2)) / (2 * variance);
////    return expf(exponent) / (2 * PI * variance);
////}
////
////inline __device__
////float gaussian1d_gpu_reg(float x, float variance, float
////sqrt_pi_variance)
////{
////    float gaussian1d = -(x*x) / (2 * variance);
////    gaussian1d = __expf(gaussian1d);
////    gaussian1d /= sqrt_pi_variance;
////    return gaussian1d;
////}
////
////float* generateGaussianKernel(int radius, float sigma)
////{
////    int area = (2 * radius + 1)*(2 * radius + 1);
////    float* res = new float[area];
////    for (int x = -radius; x <= radius; x++)
////        for (int y = -radius; y <= radius; y++)
////        {
////            //Co_to_idx inspired
////            int position = (x + radius)*(radius * 2 + 1) + y + radius;
////            res[position] = gaussian2d(x, y, sigma);
////        }
////    return res;
////}
//
//void PrintArray2D(unsigned int *h_arr, size_t img_x, size_t img_y)
//{
//   for (int k = 0; k < img_y; ++k) {
//      for (int l = 0; l < img_x; ++l) {
//         std::cout << h_arr[k*img_x + l] << " ";
//      }
//      std::cout << std::endl;
//   }
//   std::cout << std::endl;
//}
//
//// a block can hold at most 512 (below 2.0) or 1024 threads
//// shared memory within blocks
//// 16?2 so that they contain 512 threads
//// square<<<1, 64>>> = square<<<dim3(1, 1, 1), dim3(64, 1, 1), shmem in bytes>>>, dim(x, y, z)
//int main()
//{
//    const size_t img_x = 32;
//    const size_t img_y = 32;
//    //const size_t maxThreadsPerBlock = 1024;
//    const dim3 blockSize(16, 16, 1);
//   const size_t grid_x = (img_x / blockSize.x) + ((img_x % blockSize.x > 0) ? 1 : 0);
//   const size_t grid_y = (img_y / blockSize.y) + ((img_y % blockSize.y > 0) ? 1 : 0);
//   const dim3 gridSize(grid_x, grid_y, 1);
//
//    unsigned int *d_arr;
//   checkCudaErrors(cudaMalloc((void **)&d_arr, sizeof(int) * img_x * img_y));
//
//   //SetToZero << <gridSize, blockSize >> >(d_arr, img_x, img_y);
//   doThingWithCuda << <gridSize, blockSize >> >(d_arr, img_x, img_y);
//    cudaDeviceSynchronize();
//   checkCudaErrors(cudaGetLastError());
//
//   unsigned int *h_arr = new unsigned int[img_x * img_y];
//   checkCudaErrors(cudaMemcpy(h_arr, d_arr, sizeof(int) * img_x * img_y, cudaMemcpyDeviceToHost));
//   PrintArray2D(h_arr, img_x, img_y);
//
//    //const int size_x = 4;
//    //const int size_y = 4;
//    //const int size = size_x * size_y;
//    //const int a[size] = { 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44 };
//    //const int b[size] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600 };
//    //int c[size] = { 0 };
//
//    //// Add vectors in parallel.
//    //cudaError_t cudaStatus = addWithCuda(c, a, b, size_x, size_y);
//    //if (cudaStatus != cudaSuccess) {
//    //    fprintf(stderr, "addWithCuda failed!");
//    //    return 1;
//    //}
//
//    ////printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n", c[0], c[1], c[2], c[3], c[4]);
//    //for (int k = 0; k < size; ++k) {
//    //    printf("%d ", c[k]);
//    //}
//    //printf("\n");
//
//    //// cudaDeviceReset must be called before exiting in order for profiling and
//    //// tracing tools such as Nsight and Visual Profiler to show complete traces.
//    //cudaStatus = cudaDeviceReset();
//    //if (cudaStatus != cudaSuccess) {
//    //    fprintf(stderr, "cudaDeviceReset failed!");
//    //    return 1;
//    //}
//
//    return 0;
//}
