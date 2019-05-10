#ifndef _RREGION_BASE_CUH_
#define _RREGION_BASE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

namespace RegionBase
{
   class RegionBase
   {
   public:
      typedef std::unordered_set<RegionBase*> RBListT;

   protected:
      TransformableRegion * mParent;
      IMaterialScatterModelT * mScatterModel;
      ShapeT * mShape;
      RBListT mSubRegions;

   public:
      //void updateMaterial(const MaterialT& oldMat, const IMaterialScatterModelT& newMat);
      void updateMaterial(const IMaterialScatterModelT& oldMat, IMaterialScatterModelT& newMat);
      const MaterialT& getMaterial() const;
      IMaterialScatterModelT* getScatterModel() const;
      const RBListT& getSubRegions() const;
      void addRegion(RegionBase&); // protected member inaccessble via pointer
      const RegionBase* containingSubRegion(double pos[]) const;
      Element::UnorderedSetT getElements(bool recurse) const;
      const RegionBase* findEndOfStep(double p0[], double p1[]) const;
      const ShapeT* getShape() const;

   //protected:
      double getAtomsPerCubicMeter(const ElementT& el) const;
      const RegionBase* getParent() const;
      bool isContainingRegion(const RegionBase& searchTarget) const;
      char const * toString() const;
   };

   class TransformableRegion : public RegionBase, public ITransformT
   {
   public:
      void rotate(const double pivot[], double phi, double theta, double psi) override;
      void translate(const double distance[]) override;
   };

   class Region : public TransformableRegion
   {
   public:
      Region(Region* const parent, IMaterialScatterModelT * const msm, ShapeT const * const shape);
      Region(const Region&);
      void removeSubRegion(RegionBase& subRegion);
      void clearSubRegions();
   };
}

#endif