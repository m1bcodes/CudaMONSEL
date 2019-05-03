
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>

#include "Tests\HasherTest.cuh"
#include "Tests\LinkedListTest.cuh"
#include "Tests\MapTest.cuh"
#include "Tests\SetTest.cuh"

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
   if (err != cudaSuccess) {
      std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
      std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
      exit(1);
   }
}

__global__ void TestKernel()
{
   printf("INT_MIN %d\n", INT_MIN);
   printf("INT_MAX %d\n", INT_MAX);

   HasherTest::TestOne();

   LinkedListTest::LinkedListTest lltest;
   lltest.TestAddAllAsSet();

   SetTest::SetTest setTest;
   setTest.TestIntBasic();
   setTest.TestInt();
   setTest.TestInt2();
   setTest.TestString();
   setTest.TestSetOfSetOfString();

   MapTest::MapTest mapTest;
   mapTest.TestInteger();
   mapTest.TestString();
   mapTest.TestMapOfMap();
}

struct CompareDouble
{
   inline bool operator() (double& a, double& b) {
      return a < b;
   }
};

int main()
{
   //printf("%d\n", sizeof(CompareDouble));
   TestKernel<<<1, 1>>> ();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

    return 0;
}
