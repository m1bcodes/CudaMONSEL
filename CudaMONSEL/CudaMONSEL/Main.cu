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

#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

#include "gov\nist\microanalysis\EPQLibrary\EdgeEnergy.cuh"
#include "gov\nist\microanalysis\EPQLibrary\MaterialFactory.cuh"
#include "gov\nist\microanalysis\EPQLibrary\BrowningEmpiricalCrossSection.cuh"
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
   printf("%d, %d, %d, %d, %d, %d\n", threadIdx.x, blockIdx.x, threadIdx.y, blockIdx.y, threadIdx.z, blockIdx.z);
   printf("%d, %d, %d\n", gridDim.x, gridDim.y, gridDim.z);

   //int i = threadIdx.x + blockDim.x * blockIdx.x;
   Math2Test::testRandom1Cuda();
   Math2Test::testRandom2Cuda();

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

const unsigned int NUM_ROWS = 1;
const unsigned int NUM_COLS = 16;

__global__ void printRand()
{
   const unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;
   printf("thread id %d: %.10e\n", i, Math2::random());
}

void testGPU()
{
   printf("-----------------GPU-----------------------------\n");
   Math2::initCudaStates << <1, 1 >> >(NUM_ROWS * NUM_COLS);
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   printRand << <1, 16 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   testKernel << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
}

void testsCPU()
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
   VectorTest::testFour();

   StackTest::testOne();

   EdgeEnergy::DiracHartreeSlaterIonizationEnergies::loadxionUis();
   EdgeEnergy::NISTEdgeEnergy::loadNISTxrtdb();
   EdgeEnergy::ChantlerEdgeEnergy::loadFFastEdgeDB();
   EdgeEnergy::DTSAEdgeEnergy::loadEdgeEnergies();

   Element::init();
   MaterialFactory::init();

   BrowningEmpiricalCrossSection::init();
   ScreenedRutherfordScatteringAngle::init();
   CzyzewskiMottScatteringAngle::init();
   NISTMottScatteringAngle::init();
   GasScatteringCrossSection::init();
   NISTMottRS::init();
   MeanIonizationPotential::Berger64.readTabulatedValues();
   MeanIonizationPotential::Berger83.readTabulatedValues();
   //MeanIonizationPotential::Berger64MeanIonizationPotential::readTabulatedValues();
   //MeanIonizationPotential::Berger83MeanIonizationPotential::readTabulatedValues();

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

__global__ void printNullReference()
{
   printf("%s\n", Reference::dNullReference->getLongForm().c_str());
}

__device__ float *d_arr;
__device__ float d_arr1[8];

__global__ void initArr(float *d_tmp)
{
   d_arr = new float[8];
   memcpy(d_arr, d_tmp, sizeof(float) * 8);
}

__global__ void cpyker()
{
   memcpy(d_arr1, d_arr, sizeof(float) * 8);
   for (int i = 0; i < 8; ++i) {
      printf("%lf ", d_arr1[i]);
   }
   printf("\n");
}

__global__ void printSpwem()
{
   printf("GPU: %d\n", NISTMottScatteringAngle::getNISTMSA(59).getSpwem().size());
   for (auto a : NISTMottScatteringAngle::getNISTMSA(59).getSpwem()) {
      printf("%.10e ", a);
   }
   printf("GPU end\n");
}

__global__ void printMeanIonizationPotential()
{
   printf("GPU:\n");
   printf("Berger64: %d\n", MeanIonizationPotential::dBerger64->getData().size());
   for (auto a : MeanIonizationPotential::dBerger64->getData()) {
      printf("%.10e ", a);
   }
   printf("Berger64 end\n");
   printf("Berger83: %d\n", MeanIonizationPotential::dBerger83->getData().size());
   for (auto a : MeanIonizationPotential::dBerger83->getData()) {
      printf("%.10e ", a);
   }
   printf("Berger83 end\n");
   printf("GPU end\n");
}

void initCuda()
{
   printf("-----------------initCuda-----------------------------\n");

   char *d_data = nullptr;
   checkCudaErrors(cudaMalloc((void **)&d_data, sizeof(char) * 128));
   checkCudaErrors(cudaMemcpy(d_data, Reference::NullReference.getReference().c_str(), sizeof(char) * 128, cudaMemcpyHostToDevice));
   Reference::initCuda << <1, 1 >> >(d_data);
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   checkCudaErrors(cudaFree(d_data));
   printNullReference << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   Element::copyDataToDevice();
   Element::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   BrowningEmpiricalCrossSection::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   ScreenedRutherfordScatteringAngle::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   
   NISTMottScatteringAngle::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   NISTMottScatteringAngle::copyDataToCuda();
   printSpwem << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   printf("CPU: %d\n", NISTMottScatteringAngle::getNISTMSA(59).getSpwem().size());
   for (auto a : NISTMottScatteringAngle::getNISTMSA(59).getSpwem()) {
      printf("%.10e ", a);
   }
   printf("CPU end\n");
   NISTMottRS::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   MeanIonizationPotential::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   MeanIonizationPotential::copyDataToCuda();
   printMeanIonizationPotential << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   printf("CPU:\n");
   printf("Berger64: %d\n", MeanIonizationPotential::Berger64.getData().size());
   for (auto a : MeanIonizationPotential::Berger64.getData()) {
      printf("%.10e ", a);
   }
   printf("Berger64 end\n");
   printf("Berger83: %d\n", MeanIonizationPotential::Berger83.getData().size());
   for (auto a : MeanIonizationPotential::Berger83.getData()) {
      printf("%.10e ", a);
   }
   printf("Berger83 end\n");
   printf("CPU end\n");

   float tmp[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
   float *d_tmp;
   checkCudaErrors(cudaMalloc((void**)&d_tmp, sizeof(float) * 8));
   checkCudaErrors(cudaMemcpy(d_tmp, tmp, sizeof(float) * 8, cudaMemcpyHostToDevice));
   initArr << <1, 1 >> >(d_tmp);
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   cpyker << <1, 1 >> >();
   checkCudaErrors(cudaFree(d_tmp));

   //char *d_data;
   //checkCudaErrors(cudaMalloc((void **)&d_data, sizeof(char) * 128));
   //checkCudaErrors(cudaMemcpy(d_data, Reference::NullReference.getReference().c_str(), sizeof(char) * 128, cudaMemcpyHostToDevice));
   //Reference::initCuda << <1, 1 >> >(d_data);
   //checkCudaErrors(cudaDeviceSynchronize());
   //checkCudaErrors(cudaGetLastError());
   //checkCudaErrors(cudaFree(d_data));
   //print << <1, 1 >> >();
   //checkCudaErrors(cudaDeviceSynchronize());
   //checkCudaErrors(cudaGetLastError());
}

__device__ curandState state[16];

__device__ double generateRandomNumber()
{
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   const double res = curand_uniform_double(&state[id]);
   printf("thread no. %d: %.10e\n", id, res);
   return res;
}

__global__ void printRandomNumbers()
{
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   curand_init(0, id, id, &state[id]);
   printf("thread no. %d: %.10e\n", id, generateRandomNumber());
}

int main()
{
   //CudaClass c;
   //CudaClass *d_c = nullptr;
   //checkCudaErrors(cudaMalloc((void **)&d_c, sizeof(CudaClass)));
   //checkCudaErrors(cudaMemcpy(d_c, &c, sizeof(CudaClass), cudaMemcpyHostToDevice));
   ////int *hdata;
   ////checkCudaErrors(cudaMalloc((void **)&hdata, sizeof(int)));
   ////checkCudaErrors(cudaMemcpy(hdata, c.get_data(), sizeof(int), cudaMemcpyHostToDevice));
   ////checkCudaErrors(cudaMemcpy(d_c->get_data(), hdata, sizeof(int), cudaMemcpyDeviceToDevice));
   ////checkCudaErrors(cudaMemcpy(d_c->get_data(), c.get_data(), sizeof(int) * 8, cudaMemcpyHostToDevice));
   //checkCudaErrors(cudaMemcpyToSymbol(d_ptr, &d_c, sizeof(CudaClass*)));
   //useClass<<<1, 1>>>();
   //checkCudaErrors(cudaDeviceSynchronize());
   //checkCudaErrors(cudaFree(d_c));

   printRandomNumbers << <1, 16 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   testsCPU();
   testGPU();

   initCuda();

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
