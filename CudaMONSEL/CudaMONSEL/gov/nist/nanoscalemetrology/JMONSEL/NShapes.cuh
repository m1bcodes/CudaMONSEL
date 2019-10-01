#ifndef _NSHAPE_CUH_
#define _NSHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NShapes
{
   //class CustomShape
   //{
   //public:
   //   __host__ __device__ virtual void calcGroundtruth() = 0;
   //};

   struct LineParams
   {
      __host__ __device__ LineParams(const float h, const float w, const float linelength, const float thetal, const float thetar, const float radl, const float radr, const float x) :
         h(h), 
         w(w), 
         linelength(linelength), 
         thetal(thetal), 
         thetar(thetar), 
         radl(radl), 
         radr(radr), 
         x(x)
      {};
      const float h, w, linelength, thetal, thetar, radl, radr; // constructor params
      const float x; // translate
   };

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
      __host__ __device__ void calcRasterization(const PlaneT&, const double*, const double*, const float, const float, char*, const unsigned int, const unsigned int) const;

      __host__ __device__ NormalIntersectionShapeT* get();

      __host__ __device__ void addRestrainingPlanes();

   private:
      const double topz; // z of the top face
      const double width; // line width
      const double length; // length of line
      const double thetal; // angle of right sidewall
      const double thetar; // angle of left sidewall
      const double radtl; // radius of top right corner
      const double radtr; // radius of top left corner
      const double radbl; // radius of bottom right corner
      const double radbr; // radius of bottom left corner

      NormalMultiPlaneShapeT* enclosure;
      PlaneT* pl0;
      PlaneT* pl1;
      PlaneT* pl2;
      PlaneT* pl3;

      // right
      NormalShapeT* rightSide;
      PlaneT* pl4;

      // top right
      NormalMultiPlaneShapeT* toprightNMPS;
      PlaneT* pltr;
      NormalCylindricalShapeT* trcylinder;
      NormalUnionShapeT* toprightSide;

      // bottom right
      NormalMultiPlaneShapeT* botrightNMPS;
      PlaneT* plbr;
      PlaneT* plbr_0;
      PlaneT* plbr_1;
      NormalCylindricalShapeT* brcylinder;
      NormalDifferenceShapeT* botrightSide;

      // left
      NormalShapeT* leftSide;
      PlaneT* pl5;

      // top left
      NormalMultiPlaneShapeT* topleftNMPS;
      PlaneT* pltl;
      NormalCylindricalShapeT* tlcylinder;
      NormalUnionShapeT* topleftSide;

      // bot left
      NormalMultiPlaneShapeT* botleftNMPS;
      PlaneT* plbl;
      PlaneT* plbl_0;
      PlaneT* plbl_1;
      NormalCylindricalShapeT* blcylinder;
      NormalDifferenceShapeT* botleftSide;

      // together
      NormalIntersectionShapeT* nts;
      NormalIntersectionShapeT* nis; // the entire shape

      PlaneT* plrestr0;
      //PlaneT* plrestr1;
      NormalIntersectionShapeT* nis2;

      LineShapeT* gt0;
      LineShapeT* gt1;
      LineShapeT* gt2;
      LineShapeT* gt3;
   };

   class HorizontalStrip
   {
   public:
      __host__ __device__ HorizontalStrip(const double width, const bool fadeBottom = false);
      __host__ __device__ ~HorizontalStrip();

      __host__ __device__ void calcGroundtruth();
      __host__ __device__ void calcRasterization(const PlaneT&, const double*, const double*, const float, const float, char*, const unsigned int, const unsigned int) const;
      __host__ __device__ NormalMultiPlaneShapeT* get();

   private:
      NormalMultiPlaneShapeT* enclosure;
      PlaneT* bnd;
      PlaneT* pos;
      PlaneT* neg;

      LineShapeT* gt0;
      LineShapeT* gt1;
   };

   class Washer
   {
   public:
      __host__ __device__ Washer(const double innerRadius, const double outerRadius);
      __host__ __device__ ~Washer();

      __host__ __device__ void calcGroundtruth();
      __host__ __device__ void calcRasterization(const PlaneT&, const double*, const double*, const float, const float, char*, const unsigned int, const unsigned int) const;
      __host__ __device__ NormalDifferenceShapeT* get();

   private:
      NormalCylindricalShapeT* inner;
      NormalCylindricalShapeT* outer;
      NormalDifferenceShapeT* diff;
   };


   extern __host__ __device__ void TestProjection();
}

#endif