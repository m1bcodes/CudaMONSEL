/*
* - without length, the array parameter (eg const double a[]) is always 3 dimensional
* - use "auto" keyword for pointers ONLY
*/

#include <stdio.h>

#include <cuda_runtime.h>

#include "Amphibian\Tests\HasherTest.cuh"
#include "Amphibian\Tests\StringTest.cuh"
#include "Amphibian\Tests\LinkedListTest.cuh"
#include "Amphibian\Tests\SetTest.cuh"
#include "Amphibian\Tests\MapTest.cuh"
#include "Amphibian\Tests\VectorTest.cuh"
#include "Amphibian\Tests\StackTest.cuh"

#include "CudaUtil.h"
#include "ImageUtil.h"
#include <curand.h>
#include <curand_kernel.h>

#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MaterialFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MeanIonizationPotential.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\NISTMottRS.cuh"

#include "gov\nist\microanalysis\EPQTests\UncertainValue2Test.cuh"
#include "gov\nist\microanalysis\EPQTests\ElementTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MaterialTest.cuh"
#include "gov\nist\microanalysis\EPQTests\AtomicShellTest.cuh"
#include "gov\nist\microanalysis\EPQTests\EdgeEnergyTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MeanIonizationPotentialTest.cuh"
#include "gov\nist\microanalysis\EPQTests\SphereTest.cuh"
#include "gov\nist\microanalysis\EPQTests\Math2Test.cuh"
#include "gov\nist\microanalysis\EPQTests\CylindricalShapeTest.cuh"
#include "gov\nist\microanalysis\EPQTests\SumShapeTest.cuh"
#include "gov\nist\microanalysis\EPQTests\BetheElectronEnergyLossTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MonteCarloSSTest.cuh"

#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayers.cuh"

//__device__ __host__ float function(float x)
//{
//   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      return 10.0f * __sinf(x);
//   #else // host code here
//   #endif
//}

__global__ void testKernel()
{
   curandState state;
   int i = threadIdx.x + blockDim.x * blockIdx.x;
   curand_init(1234, i, 0, &state);
   Math2Test::testRandom1CUDA(state);
   Math2Test::testRandom2CUDA(state);

   HasherTest::TestOne();

   StringTest::EmptyTest();
   StringTest::TestOne();
   StringTest::AtoITest();
   StringTest::AtoFTest();
   StringTest::ItoATest();
   StringTest::findTest();
   StringTest::addTest();

   LinkedListTest::LinkedListTest lltest;
   lltest.InsertionTest();
   lltest.TestAddAllAsSet();
   LinkedListTest::TestListKV();
   LinkedListTest::testDLinkedList();

   SetTest::SetTest setTest;
   setTest.testIntBasic();
   setTest.testInt();
   setTest.testInt2();
   setTest.testString();
   //setTest.TestSetOfSetOfString();

   MapTest::MapTest mapTest;
   mapTest.testInteger();
   mapTest.testString();
   mapTest.testMapOfMap();

   VectorTest::testOne();
   VectorTest::testTwo();
   VectorTest::testThree();

   StackTest::testOne();
}

void GPUTest()
{
   printf("-----------------GPU-----------------------------\n");
   testKernel<<<1, 1>>>();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
}

void CPUTests()
{
   printf("-----------------CPU-----------------------------\n");
   HasherTest::TestOne();

   StringTest::EmptyTest();
   StringTest::TestOne();
   StringTest::AtoITest();
   StringTest::AtoFTest();
   StringTest::ItoATest();
   StringTest::findTest();
   StringTest::addTest();

   LinkedListTest::LinkedListTest lltest;
   lltest.InsertionTest();
   lltest.TestAddAllAsSet();
   LinkedListTest::TestListKV();
   LinkedListTest::testDLinkedList();

   SetTest::SetTest setTest;
   setTest.testIntBasic();
   setTest.testInt();
   setTest.testInt2();
   setTest.testString();
   //setTest.TestSetOfSetOfString();

   MapTest::MapTest mapTest;
   mapTest.testInteger();
   mapTest.testString();
   mapTest.testMapOfMap();

   VectorTest::testOne();
   VectorTest::testTwo();
   VectorTest::testThree();

   StackTest::testOne();

   EdgeEnergy::DiracHartreeSlaterIonizationEnergies::loadxionUis();
   EdgeEnergy::NISTEdgeEnergy::loadNISTxrtdb();
   EdgeEnergy::ChantlerEdgeEnergy::loadFFastEdgeDB();
   EdgeEnergy::DTSAEdgeEnergy::loadEdgeEnergies();

   Element::init();
   MaterialFactory::init();

   ScreenedRutherfordScatteringAngle::init();
   CzyzewskiMottScatteringAngle::init();
   NISTMottScatteringAngle::init();
   GasScatteringCrossSection::init();
   NISTMottRS::init();
   MeanIonizationPotential::Berger64MeanIonizationPotential::readTabulatedValues();
   MeanIonizationPotential::Berger83MeanIonizationPotential::readTabulatedValues();

   Math2Test::testRandom1();
   Math2Test::testRandom2();

   UncertainValue2::UncertainValue2 v0(0, "abc", 5);
   UncertainValue2::UncertainValue2 v1(1);
   UncertainValue2::UncertainValue2 v2(2, 10);
   UncertainValue2::UncertainValue2 v3(2, 10);
   printf("%d\n", v1.equals(v2));
   printf("%d\n", v1.equals(v3));
   printf("%d\n", v2.equals(v3));

   UncertainValue2Test::UncertainValue2Test uvTest;
   uvTest.testSpecialValues();
   uvTest.testA();
   uvTest.testB();
   uvTest.testC();
   uvTest.testAB();
   uvTest.testAdd1();
   uvTest.testAdd2();
   uvTest.testAdd3();
   uvTest.testMultiply();
   uvTest.testDivide();
   uvTest.testFunctions();

   ElementTest::ElementTest elementTest;
   elementTest.testZero();
   elementTest.testOne();

   MaterialTest::MaterialTest mat;
   mat.testOne();

   AtomicShellTest::testOne();

   EdgeEnergyTest::testOne();

   MeanIonizationPotentialTest::testOne();

   SphereTest::testContains();
   SphereTest::testGetFirstIntersection();

   CylindricalShapeTest::CylindricalShapeTest cylindricalShapeTest;
   cylindricalShapeTest.testZero();
   cylindricalShapeTest.testOne();
   cylindricalShapeTest.testTwo();
   cylindricalShapeTest.testThree();
   cylindricalShapeTest.testFour();
   cylindricalShapeTest.testFive();
   cylindricalShapeTest.testSix();
   cylindricalShapeTest.testSeven();
   cylindricalShapeTest.testEight();
   cylindricalShapeTest.testNine();
   cylindricalShapeTest.testTen();
   cylindricalShapeTest.testEleven();
   cylindricalShapeTest.testTwelve();

   BetheElectronEnergyLossTest::testOne();

   MonteCarloSSTest::testOne();

   SumShapeTest::SumShapeTest sumShapeTest;
   sumShapeTest.testGetFirstIntersection();
   sumShapeTest.testAll();
}

__device__ int n = 0;

__global__ void printKernel(const float *a, int numElements)
{
   for (int i = 0; i < numElements; ++i) {
      printf("%lf ", a[i]);
   }
   printf("done\n");

   n += 1;
   printf("%d\n", n);

   __syncthreads();
}

__global__ void vectorAdd(const float *A, const float *B, float *C, int numElements)
{
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   if (i < numElements)
   {
      C[i] = A[i] + B[i];
   }
}

__global__ void wewKernel()
{
   //Reference::dNullReference = new Reference::CrudeReference("WEW");
   //memcpy(Reference::dNullReference, tmp, sizeof(Reference::CrudeReference));
   printf("123\n");
   //printf("%s\n", tmp->getLongForm().c_str());
}

int main()
{
   int numElements = 8;
   size_t size = numElements * sizeof(float);
   printf("[Vector addition of %d elements]\n", numElements);

   float *h_A = (float *)malloc(size);
   float *h_B = (float *)malloc(size);
   float *h_C = (float *)malloc(size);

   if (h_A == NULL || h_B == NULL || h_C == NULL) {
      fprintf(stderr, "Failed to allocate host vectors!\n");
      exit(EXIT_FAILURE);
   }

   for (int i = 0; i < numElements; ++i) {
      h_A[i] = i;
      h_B[i] = i;
   }

   float *d_A = NULL;
   checkCudaErrors(cudaMalloc((void **)&d_A, size));

   float *d_B = NULL;
   checkCudaErrors(cudaMalloc((void **)&d_B, size));

   float *d_C = NULL;
   checkCudaErrors(cudaMalloc((void **)&d_C, size));

   printKernel << <1, 1 >> >(d_A, 8);
   checkCudaErrors(cudaGetLastError());
   printKernel << <1, 1 >> >(d_B, 8);
   checkCudaErrors(cudaGetLastError());
   printKernel << <1, 1 >> >(d_C, 8);
   checkCudaErrors(cudaGetLastError());

   //checkCudaErrors(cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice));
   //checkCudaErrors(cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice));

   //int threadsPerBlock = 4;
   //int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
   //vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, numElements);
   //checkCudaErrors(cudaGetLastError());
   //checkCudaErrors(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

   //for (int i = 0; i < numElements; ++i) {
   //   if (fabs(h_A[i] + h_B[i] - h_C[i]) > 1e-5) {
   //      fprintf(stderr, "Result verification failed at element %d!\n", i);
   //      exit(EXIT_FAILURE);
   //   }
   //}
   //printf("Test PASSED\n");

   checkCudaErrors(cudaFree(d_A));
   checkCudaErrors(cudaFree(d_B));
   checkCudaErrors(cudaFree(d_C));

   // Free host memory
   free(h_A);
   free(h_B);
   free(h_C);

   printf("Done\n");

   //Reference::CrudeReference hcr("ABCDEFG"), hcr1("");
   //Reference::CrudeReference *tmp = nullptr;
   //checkCudaErrors(cudaMalloc((void **)&tmp, sizeof(Reference::CrudeReference)));
   //checkCudaErrors(cudaMemcpy(tmp, &hcr, sizeof(Reference::CrudeReference), cudaMemcpyHostToDevice));
   //checkCudaErrors(cudaMemcpyToSymbol(Reference::dNullReference, tmp, sizeof(Reference::CrudeReference)));
   wewKernel<<<1, 1>>>();
   checkCudaErrors(cudaGetLastError());
   //checkCudaErrors(cudaMemcpy(&hcr1, Reference::dNullReference, sizeof(Reference::CrudeReference), cudaMemcpyDeviceToHost));
   //checkCudaErrors(cudaFree(tmp));
   //printf("%s\n", hcr1.getLongForm().c_str());
   //tmp = (Reference::CrudeReference*)malloc(sizeof(Reference::CrudeReference));
   //checkCudaErrors(cudaMemcpy(tmp, &hcr, sizeof(Reference::CrudeReference), cudaMemcpyDeviceToHost));
   //checkCudaErrors(cudaMemcpyFromSymbol((void**)&tmp, "dNullReference", sizeof(Reference::CrudeReference), 0, cudaMemcpyDeviceToHost));
   //printf("%s\n", tmp->getLongForm().c_str());
   //free(tmp);

   CPUTests();
   GPUTest();

   LinesOnLayers::run();

   return 0;
}

//static const int N = 2;
//
//class vecarray
//{
//public:
//   int *vecptr[N]; //array of pointers pointing to array
//   int dim[N]; //store length of each array pointed to
//
//   __device__ __host__ vecarray(); //constructor
//   __device__ __host__ int sum();  //sum up all the elements in the array being pointed to
//};
//
//vecarray::vecarray()
//{
//   for (int i = 0; i<N; i++)
//   {
//      vecptr[i] = NULL;
//      dim[i] = 0;
//   }
//}
//
//__device__ __host__ int vecarray::sum()
//{
//   int i = 0, j = 0, s = 0;
//   for (i = 0; i<N; i++)
//      for (j = 0; j < dim[i]; j++)
//         s += vecptr[i][j];
//   return s;
//}
//
//__global__ void addvecarray(vecarray * v, int *s)
//{
//   *s = v->sum();
//}
//
//int main()
//{
//   //copy *V to device, do sum() and pass back
//   vecarray *v, *dev_v; // the result by dev_v
//   v = new vecarray;
//   int a[3] = { 1, 2, 3 }; //initialize v manually
//   int b[4] = { 4, 5, 6, 7 };
//   int result = 0;
//   int *dev_result;
//   v->vecptr[0] = a;
//   v->vecptr[1] = b;
//   v->dim[0] = 3; v->dim[1] = 4;
//   int *vptr[N];
//
//   checkCudaErrors(cudaMalloc((void**)&dev_v, sizeof(vecarray)));
//   checkCudaErrors(cudaMemcpy(dev_v, v, sizeof(vecarray), cudaMemcpyHostToDevice)); //copy class object
//   checkCudaErrors("cudaMemcpy1 fail");
//
//   for (int i = 0; i < N; i++){
//      checkCudaErrors(cudaMalloc((void**)&(vptr[i]), v->dim[i] * sizeof(int)));
//      checkCudaErrors(cudaMemcpy(&(dev_v->vecptr[i]), &vptr[i], sizeof(int*), cudaMemcpyHostToDevice));
//   }
//
//   for (int i = 0; i<N; i++) { // copy arrays
//      checkCudaErrors(cudaMemcpy(vptr[i], v->vecptr[i], v->dim[i] * sizeof(int), cudaMemcpyHostToDevice));
//   }
//   checkCudaErrors(cudaMalloc((void **)&dev_result, sizeof(int)));
//   addvecarray<<<1, 1>>>(dev_v, dev_result);
//
//   checkCudaErrors(cudaMemcpy(&result, dev_result, sizeof(int), cudaMemcpyDeviceToHost));
//   printf("the result is %d\n", result);
//
//   return 0;
//}
