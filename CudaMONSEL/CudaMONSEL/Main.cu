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

#include "gov\nist\microanalysis\EPQLibrary\CzyzewskiMottScatteringAngle.cuh"
#include "gov\nist\microanalysis\EPQLibrary\GasScatteringCrossSection.cuh"
#include "gov\nist\microanalysis\EPQLibrary\NISTMottScatteringAngle.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\ScreenedRutherfordScatteringAngle.cuh"
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

#include "ImageUtil.h"

int main()
{
   CzyzewskiMottScatteringAngle::init();
   NISTMottScatteringAngle::init();
   GasScatteringCrossSection::init();
   NISTMottRS::init();

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

   CylindricalShapeTest::testZero();
   CylindricalShapeTest::testOne();
   CylindricalShapeTest::testTwo();
   CylindricalShapeTest::testThree();
   CylindricalShapeTest::testFour();
   CylindricalShapeTest::testFive();
   CylindricalShapeTest::testSix();
   CylindricalShapeTest::testSeven();
   CylindricalShapeTest::testEight();
   CylindricalShapeTest::testNine();
   CylindricalShapeTest::testTen();
   CylindricalShapeTest::testEleven();
   CylindricalShapeTest::testTwelve();

   //BetheElectronEnergyLossTest::testOne();

   //SumShapeTest::testGetFirstIntersection();
   //SumShapeTest::testAll();

   return 0;
}
