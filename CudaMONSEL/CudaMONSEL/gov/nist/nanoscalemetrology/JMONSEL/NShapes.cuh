#ifndef _NSHAPE_CUH_
#define _NSHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NShapes
{
   struct Line
   {
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
      NormalIntersectionShapeT* nis;
   };

   __host__ __device__ extern NormalShapeT* createLine(
      Line&,
      double topz, // z of the top face
      double width, // line width
      double length, // length of line
      double thetal, // angle of right sidewall
      double thetar, // angle of left sidewall
      double radl, // radius of top right corner
      double radr // radius of top left corner
      );

   __host__ __device__ extern void destroyLine(Line&);
}

#endif