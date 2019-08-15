#ifndef _MULTI_PLANE_SHAPE_CUH_
#define _MULTI_PLANE_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQlibrary\ITransform.cuh"

#include "Amphibian\vector.cuh"

namespace MultiPlaneShape
{
   class Plane : public ShapeT, public ITransformT
   {
   public:
      __host__ __device__ Plane(const double normal[], const double point[]);
      Plane(const Plane&);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      bool almostContains(const double p[]);

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      __host__ __device__ const double* getNormal() const;
      __host__ __device__ const double* getPoint() const;

   private:
      double mNormal[3];
      double mPoint[3];
   };

   class MultiPlaneShape : public ShapeT, public ITransformT//, TrajectoryVRML.IRender
   {
   public:
      typedef amp::vector<Plane*> Planes;

      __host__ __device__ MultiPlaneShape();
      __host__ __device__ MultiPlaneShape(Plane* const planes[], int len);

      __host__ __device__ bool contains(const double pos[]) const override;
      __host__ __device__ double getFirstIntersection(const double pos0[], const double pos1[]) override;
      __host__ __device__ StringT toString() const override;

      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;

      //void addOffsetPlane(const double normal[], const double pt[], double dist);
      //void addPlane(const double normal[], const double point[]);
      __host__ __device__ virtual void addPlane(Plane *plane);

      Planes getPlanes() const;

   protected:
      Planes mPlanes;
      double mInsidePos[3];
   };

   //extern MultiPlaneShape createFilm(const double normal[], const double pt1[], double thickness);
   //extern MultiPlaneShape createSubstrate(const double normal[], const double pt[]);
   //extern MultiPlaneShape createBlock(const double dims[], const double point[], double phi, double theta, double psi);
   //extern MultiPlaneShape createNamid(double center[], int n, double height, double base);

   struct LineShape
   {
      __host__ __device__ LineShape();
      __host__ __device__ LineShape(const double*, const double*);

      double P0[3];
      double P1[3];
   };

   extern __host__ __device__ int intersect3D_SegmentPlane(const LineShape& S, const Plane& Pn, double* I);
   extern __host__ __device__ int intersect3D_2Planes(const Plane& Pn1, const Plane& Pn2, LineShape& L);
}

#endif
