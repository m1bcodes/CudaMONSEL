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
#include "Amphibian\random.cuh"
#include "Amphibian\Algorithm.cuh"

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
#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"

#include "gov\nist\microanalysis\EPQTests\UncertainValue2Test.cuh"
#include "gov\nist\microanalysis\EPQTests\ElementTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MaterialTest.cuh"
#include "gov\nist\microanalysis\EPQTests\AtomicShellTest.cuh"
#include "gov\nist\microanalysis\EPQTests\EdgeEnergyTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MeanIonizationPotentialTest.cuh"
#include "gov\nist\microanalysis\EPQTests\SphereTest.cuh"
#include "gov\nist\microanalysis\EPQTests\CylindricalShapeTest.cuh"
#include "gov\nist\microanalysis\EPQTests\SumShapeTest.cuh"
#include "gov\nist\microanalysis\EPQTests\BetheElectronEnergyLossTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MonteCarloSSTest.cuh"

#include "gov\nist\nanoscalemetrology\JMONSELTests\LinesOnLayers0.cuh"

#include <chrono>

//__device__ __host__ float function(float x)
//{
//   #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//      return 10.0f * __sinf(x);
//   #else // host code here
//   #endif
//}

__host__ __device__ void testRandom1()
{
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());

   printf("Math2Test::testRandom1() completed.\n");
}

__host__ __device__ void testRandom2()
{
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());
   printf("%.10e\n", Random::random());

   printf("Math2Test::testRandom2() completed.\n");
}

__global__ void testLibraryCuda()
{
   printf("%d, %d, %d, %d, %d, %d\n", threadIdx.x, blockIdx.x, threadIdx.y, blockIdx.y, threadIdx.z, blockIdx.z);
   printf("%d, %d, %d\n", gridDim.x, gridDim.y, gridDim.z);

   testRandom1();
   testRandom2();

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

//const unsigned int H = 128, W = 80;
//const unsigned int TX = 16, TY = 16;
const unsigned int H = 4, W = 4;
const unsigned int TX = 2, TY = 2;
dim3 blockSize(TX, TY); // Equivalent to dim3 blockSize(TX, TY, 1);
int bx = (W + blockSize.x - 1) / blockSize.x;
int by = (H + blockSize.y - 1) / blockSize.y;
dim3 gridSize = dim3(bx, by);

__global__ void printRand()
{
   const unsigned int c = threadIdx.x + blockDim.x * blockIdx.x;
   const unsigned int r = blockIdx.y*blockDim.y + threadIdx.y;

   int blockId = blockIdx.x + blockIdx.y * gridDim.x;
   int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
   printf("%d, %d (%d): %.5e | ", r, c, threadId, Random::random());
}

void testGPU()
{
   printf("-----------------GPU-----------------------------\n");
   Random::initCudaStates << <1, 1 >> >(H * W);
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   printRand << <gridSize, blockSize >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   testLibraryCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
}

void testLibrary()
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

   testRandom1();
   testRandom2();
}

void initSim()
{
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
}

void testSim()
{
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

   //SumShapeTest::SumShapeTest sumShapeTest;
   //sumShapeTest.testGetFirstIntersection();
   //sumShapeTest.testAll();
}

__global__ void printNullReference()
{
   printf("%s\n", Reference::d_NullReference->getLongForm().c_str());
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
   printf("Berger64: %d\n", MeanIonizationPotential::d_Berger64->getData().size());
   for (auto a : MeanIonizationPotential::d_Berger64->getData()) {
      printf("%.10e ", a);
   }
   printf("Berger64 end\n");
   printf("Berger83: %d\n", MeanIonizationPotential::d_Berger83->getData().size());
   for (auto a : MeanIonizationPotential::d_Berger83->getData()) {
      printf("%.10e ", a);
   }
   printf("Berger83 end\n");
   printf("GPU end\n");
}

void initCuda()
{
   printf("-----------------initCuda-----------------------------\n");

   Material::initCuda << <1, 1 >> >();
   AlgorithmUser::initCuda << <1, 1 >> >();
   NISTMottRS::initFactory << <1, 1 >> >();
   NUTableInterpolation::initFactory << <1, 1 >> >();

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

   Element::copyDataToCuda();
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
   NISTMottScatteringAngle::initFactory << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   NISTMottScatteringAngle::transferDataToCuda();
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
   MeanIonizationPotential::transferDataToCuda();
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
}
// causes stack overflow on GPU since it is sorting on worst case (ie already sorted array)
// use case: Histogram::Histogram
//__global__ void AlgorithmTest()
//{
//   VectorXd bins(360 + 1, 0);
//   bins[0] = 0;
//   const double delta = (Math2::PI * 2) / (360 + 1);
//   for (int i = 1; i < bins.size(); ++i)
//      bins[i] = i * delta;
//   Algorithm::quicksort(bins.data(), 0, bins.size() - 1);
//}

__global__ void printTotalCrossSection()
{
   printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dC).totalCrossSection(8.01088e-17));
   //printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dO));
   //printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dH));

   //printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dC).totalCrossSection(8.01088e-17));
   //printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dO).totalCrossSection(8.01088e-17));
   //printf("%.5e\n", NISTMottScatteringAngle::d_Factory->get(*Element::dH).totalCrossSection(8.01088e-17));
}

int main()
{
   cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1e9);
   cudaDeviceSetLimit(cudaLimitStackSize, 131072);
   size_t pValue;
   cudaDeviceGetLimit(&pValue, cudaLimitMallocHeapSize);
   printf("cudaLimitMallocHeapSize: %d\n", pValue);
   cudaDeviceGetLimit(&pValue, cudaLimitStackSize);
   printf("cudaLimitStackSize: %d\n", pValue);

   testLibrary();
   initSim();
   testSim();
   testGPU();

   initCuda();

   LinesOnLayers::transferDataToCuda();

   size_t a, t;
   checkCudaErrors(cudaMemGetInfo(&a, &t));
   printf("free/total: %lu/%lu\n", a, t);

   LinesOnLayers::initCuda << <1, 1 >> >();
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());

   //LinesOnLayers::runCuda << <1, 1 >> >();
   //checkCudaErrors(cudaDeviceSynchronize());
   //checkCudaErrors(cudaGetLastError());

   //cudaStream_t streams[2];
   //cudaStreamCreate(&streams[0]);
   //cudaStreamCreate(&streams[1]);

   auto start = std::chrono::system_clock::now();
   //LinesOnLayers::runCudaSinglePixel << <gridSize, blockSize >> >();
   //LinesOnLayers::runCudaSinglePixel << <1, 1, 0, streams[0]>> >(0, 0);
   //LinesOnLayers::runCudaSinglePixel << <1, 1, 0, streams[1] >> >(0, 1);
   //LinesOnLayers::runCudaSinglePixel << <8, 1 >> >();
   LinesOnLayers::runCudaSinglePixel << <1, 1 >> >(0, 4);
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError());
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);
   std::cout << std::endl << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

   printf("done\n");

   //testsCPU();
   //LinesOnLayers::initCuda();
   //LinesOnLayers::runCuda();

   return 0;
}
