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
      NormalMultiPlaneShape();

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      bool contains(const double pos0[], const double pos1[]) const  override;
      VectorXd getFirstNormal(const double pos0[], const double pos1[]) override;
      VectorXd getPreviousNormal() const override;

      //void addPlane(const double normal[], const double point[]);
      void addPlane(PlaneT&) override;

      int getNumPlanes() const;
      VectorXd getNormal(int index) const;
      double getB(int index) const;

   private:
      void updateCach();

      VectorXd result;

      MultiPlaneShape::Planes mPlanes;

      int numPlanes; // # of planes in this shape

      /*
      * narray[i] is the normal vector of the ith plane
      */
      MatrixXd narray;

      /*
      * carray[i] is the point contained in the ith plane
      */
      MatrixXd carray;
   };

   // extern void addPlane(const double p1[], const double p2[], const double p3[]);
}

#endif