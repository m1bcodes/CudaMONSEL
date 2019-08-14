#ifndef _NSHAPE_CUH_
#define _NSHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NShapes
{
   class Line
   {
   public:
      __host__ __device__ Line(
         double topz,
         double width,
         double length,
         double thetal,
         double thetar,
         double radl,
         double radr
         );
      __host__ __device__ ~Line();

      __host__ __device__ void create();
      __host__ __device__ void calcGroundtruth();
      __host__ __device__ void calcRasterization(const PlaneT&, const double*, const double*, const float, const float);

      __host__ __device__ NormalIntersectionShapeT* get();

   private:
      double topz; // z of the top face
      double width; // line width
      double length; // length of line
      double thetal; // angle of right sidewall
      double thetar; // angle of left sidewall
      double radl; // radius of top right corner
      double radr; // radius of top left corner

      NormalMultiPlaneShapeT* enclosure;
      PlaneT* pl0;
      PlaneT* pl1;
      PlaneT* pl2;
      PlaneT* pl3;

      NormalMultiPlaneShapeT* rightNMPS;
      NormalShapeT* rightSide;
      PlaneT* pl4;
      PlaneT* plr;
      NormalCylindricalShapeT* rcylinder;

      NormalMultiPlaneShapeT* leftNMPS;
      NormalShapeT* leftSide;
      PlaneT* pl5;
      PlaneT* pll;
      NormalCylindricalShapeT* lcylinder;

      NormalIntersectionShapeT* nts;
      NormalIntersectionShapeT* nis; // the entire shape

      LineShapeT* segment0;
      LineShapeT* segment1;
      LineShapeT* segment2;
      LineShapeT* segment3;
   };

   extern __host__ __device__ void TestProjection();
}

#endif