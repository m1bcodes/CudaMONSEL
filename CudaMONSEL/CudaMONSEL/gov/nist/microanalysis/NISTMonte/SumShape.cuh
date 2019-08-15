// file: gov\nist\microanalysis\NISTMonte\SumShape.cuh

#ifndef _SUM_SHAPE_CUH_
#define _SUM_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

#include "Amphibian\vector.cuh"

namespace SumShape
{
   class SumShape : public ShapeT, public ITransformT//, TrajectoryVRML.IRender
   {
   public:
      SumShape();
      SumShape(ShapeT* const shapes[], int len);
      __host__ __device__ SumShape(ShapeT* const a, ShapeT* const b);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ const amp::vector<ShapeT*>& getShapes() const;

      void addShape(ShapeT* b);

   private:
      amp::vector<ShapeT*> mShapes;
   };

   extern __host__ __device__ void rotateShape(const double pivot[], double phi, double theta, double psi, ShapeT* shape);
   extern __host__ __device__ void translateShape(const double distance[], ShapeT* shape);
}

#endif