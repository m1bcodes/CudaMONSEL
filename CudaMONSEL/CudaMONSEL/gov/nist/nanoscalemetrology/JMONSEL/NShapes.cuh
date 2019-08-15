#ifndef _NSHAPE_CUH_
#define _NSHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NShapes
{
   class Line
   {
   public:
      __host__ __device__ Line(
         const double topz,
         const double width,
         const double length,
         const double thetal,
         const double thetar,
         const double radl,
         const double radr
         );
      __host__ __device__ ~Line();

      __host__ __device__ void create();
      __host__ __device__ void calcGroundtruth();
      __host__ __device__ void calcRasterization(const PlaneT&, const double*, const double*, const float, const float, char*, const unsigned int, const unsigned int);

      __host__ __device__ NormalIntersectionShapeT* get();

   private:
      const double topz; // z of the top face
      const double width; // line width
      const double length; // length of line
      const double thetal; // angle of right sidewall
      const double thetar; // angle of left sidewall
      const double radl; // radius of top right corner
      const double radr; // radius of top left corner

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