//#ifndef _RREGION_BASE_CUH_
//#define _RREGION_BASE_CUH_
//
//#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
//#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
////#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
//#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
//
//class RegionBase
//{
//public:
//   typedef std::vector<RegionBase*> RBListT;
//
//protected:
//   TransformableRegion* mParent;
//   IMaterialScatterModel* mScatterModel;
//   Shape* mShape;
//   RBListT mSubRegions;
//
//public:
//   RegionBase(const RegionBase& rb);
//   bool operator==(const RegionBase& rb);
//   void updateMaterial(Material::Material& oldMat, IMaterialScatterModel& newMat);
//   void updateMaterial(IMaterialScatterModel& oldMat, IMaterialScatterModel& newMat);
//   Material::Material& getMaterial();
//   IMaterialScatterModel& getScatterModel();
//   RBListT getSubRegions();
//
//protected:
//   RegionBase* containingSubRegion(double pos[]);
//   Element::UnorderedSetT getElements(bool recurse);
//   RegionBase* findEndOfStep(double p0[], double p1[]);
//   Shape& getShape();
//   double getAtomsPerCubicMeter(Element::Element& el);
//   RegionBase& getParent();
//   bool isContainingRegion(RegionBase& searchTarget);
//   char const * const toString();
//};
//
//class TransformableRegion : public RegionBase, public ITransform
//{
//public:
//   void rotate(double pivot[], double phi, double theta, double psi) override;
//   void translate(double distance[]) override;
//};
//
//#endif