#ifndef _MULTI_PLANE_SHAPE_CUH_
#define _MULTI_PLANE_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQlibrary\ITransform.cuh"

namespace MultiPlaneShape
{
   class Plane : public ShapeT, public ITransformT
   {
   public:
      Plane(const double normal[], int normallen, const double point[], int pointlen);
      Plane(const Plane&);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      bool almostContains(double p[]);

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      VectorXd getNormal() const;
      VectorXd getPoint() const;

      const VectorXd& getNormalRef() const;
      const VectorXd& getPointRef() const;

   private:
      VectorXd mNormal;
      VectorXd mPoint;
   };

   class MultiPlaneShape : public ShapeT, public ITransformT//, TrajectoryVRML.IRender
   {
   public:
      typedef std::vector<Plane*> Planes;

      MultiPlaneShape();
      MultiPlaneShape(Plane* const planes[], int len);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;

      //void addOffsetPlane(const double normal[], const double pt[], double dist);
      //void addPlane(const double normal[], const double point[]);
      virtual void addPlane(Plane& plane);

      Planes getPlanes() const;

   protected:
      Planes mPlanes;
      VectorXd mInsidePos;
   };

   //extern MultiPlaneShape createFilm(const double normal[], const double pt1[], double thickness);
   //extern MultiPlaneShape createSubstrate(const double normal[], const double pt[]);
   //extern MultiPlaneShape createBlock(const double dims[], const double point[], double phi, double theta, double psi);
   //extern MultiPlaneShape createNamid(double center[], int n, double height, double base);
}

#endif
