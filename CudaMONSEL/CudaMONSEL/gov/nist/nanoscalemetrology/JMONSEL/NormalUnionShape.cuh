//package gov.nist.nanoscalemetrology.JMONSEL;

#ifndef _NORMAL_UNION_SHAPE_CUH_
#define _NORMAL_UNION_SHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace NormalUnionShape
{
   class NormalUnionShape : public SumShapeT, public NormalShapeT
   {
   public:
      NormalUnionShape(NormalShapeT& a, NormalShapeT& b);

      bool contains(const double pos[]) const override;
      double getFirstIntersection(const double pos0[], const double pos1[]) override;
      StringT toString() const override;

      bool contains(const double pos0[], const double pos1[]) const override;
      const double* getFirstNormal(const double pos0[], const double pos1[]) override;
      const double* getPreviousNormal() const override;

      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;
      
   private:
      double result[4];
   };
}

#endif