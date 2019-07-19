// file: gov\nist\microanalysis\EPQTests\SumShapeTest.cuh

#ifndef _SUM_SHAPE_TEST_CUH_
#define _SUM_SHAPE_TEST_CUH_

#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\BasicMaterialModel.cuh"

namespace SumShapeTest
{
   class SumShapeTest
   {
   public:
      SumShapeTest();

      void testGetFirstIntersection();
      void testAll();

   private:
      void pointInside(double res[]);

      MaterialT mat1, mat2;
      SphereT sphere;
      NullMaterialScatterModelT NULL_MSM;
      RegionT mChamber;
      GaussianBeamT beam;
      MonteCarloSST mMonte;

      PlaneT pl;
      //PlaneT* pls[1];
      MultiPlaneShapeT blk;

      SphereT mCap0Outer, mCap1Outer;
      CylindricalShapeT mCylOuter;

      //ShapeT* const shapes[3];
      SumShapeT mPillOuter;

      BasicMaterialModelT bmm1, bmm2;
      RegionT r1, r2;
   };
}

#endif