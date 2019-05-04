#include <stdio.h>

#include <cuda_runtime.h>

//#include "gov/nist/microanalysis/NISTMonte/MonteCarloSS.cu"
////#include "gov/nist/microanalysis/NISTMonte/Electron.cu"
////#include "gov/nist/microanalysis/Utility/CSVReader.h"
//#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\Composition.cuh"
//#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"

#include "gov\nist\microanalysis\Utility\UncertainValue2.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

#include "gov\nist\microanalysis\EPQTests\UncertainValue2Test.cuh"
#include "gov\nist\microanalysis\EPQTests\ElementTest.cuh"
#include "gov\nist\microanalysis\EPQTests\MaterialTest.cuh"
#include "gov\nist\microanalysis\EPQTests\AtomicShellTest.cuh"

#include "ImageUtil.h"

int main()
{
   Element::InitializeElements();

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

   return 0;
}
