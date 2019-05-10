//#ifndef _NORMAL_MULTIPLANE_SHAPE_CUH_
//#define _NORMAL_MULTIPLANE_SHAPE_CUH_
//
//#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
//#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"
//
//namespace NormalMultiPlaneShape
//{
//   class NormalMultiPlaneShap : public MultiPlaneShapeT, public NormalShapeT
//   {
//   public:
//      NormalMultiPlaneShap();
//
//      bool contains(const double pos[]) const override;
//      double getFirstIntersection(const double pos0[], const double pos1[]) override;
//      char const * toString() const override;
//
//      void rotate(const double pivot[], double phi, double theta, double psi) override;
//      void translate(const double distance[]) override;
//
//      bool contains(const double pos0[], const double pos1[]) const  override;
//      VectorXd getFirstNormal(const double pos0[], const double pos1[]) const  override;
//      VectorXd getPreviousNormal() const  override;
//   };
//}
//
//#endif