#ifndef _NORMAL_MULTIPLANE_SHAPE_CUH_
#define _NORMAL_MULTIPLANE_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace NormalMultiPlaneShape
{
   class NormalMultiPlaneShape : public MultiPlaneShapeT, public NormalShapeT
   {
   public:
      __host__ __device__ NormalMultiPlaneShape();

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ bool contains(const double pos0[], const double pos1[]) const  override;
      __host__ __device__ const double* getFirstNormal(const double pos0[], const double pos1[]) override;
      __host__ __device__ const double* getPreviousNormal() const override;

      //void addPlane(const double normal[], const double point[]);
      __host__ __device__ void addPlane(PlaneT*) override;

      int getNumPlanes() const;
      //const VectorXd& getNormal(int index) const;
      double getB(int index) const;

   private:
      __host__ __device__ void updateCach();

      double result[4];

      MatrixXd narray;
      MatrixXd carray;
   };
}

#endif