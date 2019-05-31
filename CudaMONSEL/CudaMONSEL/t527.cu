//#include <stdio.h>
//#include <curand.h>
//#include <curand_kernel.h>
//#include <math.h>
//#include <assert.h>
//
//#include "CudaUtil.h"
//
//#define MIN 2
//#define MAX 7
//#define ITER 10000000
//
////__global__ void generate_kernel(curandState *my_curandstate, const unsigned int n, const unsigned *min_rand_int, const unsigned *max_rand_int, unsigned int *result)
////{
////   int idx = threadIdx.x + blockDim.x * blockIdx.x;
////   curand_init(1234, idx, 0, &my_curandstate[idx]);
////
////   for (int c = 0; c < n; ++c) {
////      int myrand = (int)truncf(randfloat(my_curandstate + idx) * (max_rand_int[idx] - min_rand_int[idx] + 0.999999) + min_rand_int[idx]);
////
////      assert(myrand <= max_rand_int[idx]);
////      assert(myrand >= min_rand_int[idx]);
////      result[myrand - min_rand_int[idx]]++;
////   }
////}
//
//__global__ void generate_kernel(const unsigned int n, const unsigned *min_rand_int, const unsigned *max_rand_int, unsigned int *result)
//{
//   curandState state;
//   int i = threadIdx.x + blockDim.x * blockIdx.x;
//   curand_init(1234, i, 0, &state);
//
//   for (int c = 0; c < n; ++c) {
//      int myrand = (int)truncf(curand_uniform(&state)* (max_rand_int[0] - min_rand_int[0] + 1) + min_rand_int[0]);
//
//      assert(myrand <= max_rand_int[0]);
//      assert(myrand >= min_rand_int[0]);
//      result[myrand - min_rand_int[0] + i * (max_rand_int[0] - min_rand_int[0] + 1)]++;
//   }
//}
//
//int main()
//{
//   //curandState *d_state;
//   //checkCudaErrors(cudaMalloc(&d_state, sizeof(curandState)));
//   unsigned int *d_result, *h_result;
//   unsigned int *d_max_rand_int, *d_min_rand_int;
//   checkCudaErrors(cudaMalloc(&d_result, 2 * (MAX - MIN + 1) * sizeof(unsigned int)));
//   h_result = (unsigned int *)malloc(2 * (MAX - MIN + 1) * sizeof(unsigned int));
//   checkCudaErrors(cudaMalloc(&d_max_rand_int, sizeof(unsigned int)));
//   checkCudaErrors(cudaMalloc(&d_min_rand_int, sizeof(unsigned int)));
//
//   checkCudaErrors(cudaMemset(d_result, 0, 2 * (MAX - MIN + 1) * sizeof(unsigned int)));
//
//   unsigned int h_min_rand_int = MIN;
//   unsigned int h_max_rand_int = MAX;
//   checkCudaErrors(cudaMemcpy(d_min_rand_int, &h_min_rand_int, sizeof(unsigned int), cudaMemcpyHostToDevice));
//   checkCudaErrors(cudaMemcpy(d_max_rand_int, &h_max_rand_int, sizeof(unsigned int), cudaMemcpyHostToDevice));
//   //generate_kernel<<<1, 1>>>(d_state, ITER, d_min_rand_int, d_max_rand_int, d_result);
//   generate_kernel<<<1, 2>>>(ITER, d_min_rand_int, d_max_rand_int, d_result);
//   checkCudaErrors(cudaMemcpy(h_result, d_result, 2 * (MAX - MIN + 1) * sizeof(unsigned int), cudaMemcpyDeviceToHost));
//   printf("Bin:    Count: \n");
//   for (int i = 0; i < 2 * (MAX - MIN + 1); i++)
//      printf("%d    %d\n", i, h_result[i]);
//
//   checkCudaErrors(cudaFree(d_min_rand_int));
//   checkCudaErrors(cudaFree(d_max_rand_int));
//   free(h_result);
//   checkCudaErrors(cudaFree(d_result));
//
//   return 0;
//}