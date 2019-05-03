#ifndef _RREGION_BASE_CUH_
#define _RREGION_BASE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

#include <vector>

namespace RegionBase
{
   class RegionBase
   {
   public:
      typedef std::vector<RegionBase*> RBListT;

   protected:
      TransformableRegion * mParent;
      IMaterialScatterModelT const * mScatterModel;
      ShapeT const * mShape;
      RBListT mSubRegions;

   public:
      RegionBase();
      RegionBase(const RegionBase& rb);
      bool operator==(const RegionBase& rb) const;
      void updateMaterial(const MaterialT& oldMat, const IMaterialScatterModelT& newMat);
      void updateMaterial(const IMaterialScatterModelT& oldMat, const IMaterialScatterModelT& newMat);
      const MaterialT& getMaterial() const;
      const IMaterialScatterModelT& getScatterModel() const;
      RBListT getSubRegions() const;
      void addRegion(RegionBase * const);

   protected:
      const RegionBase* containingSubRegion(double pos[]) const;
      Element::UnorderedSetT getElements(bool recurse) const;
      const RegionBase* findEndOfStep(double p0[], double p1[]) const;
      const ShapeT& getShape() const;
      double getAtomsPerCubicMeter(const ElementT& el) const;
      const RegionBase& getParent() const;
      bool isContainingRegion(const RegionBase& searchTarget) const;
      char const * toString() const;
   };

   class TransformableRegion : public RegionBase, public ITransformT
   {
   public:
      void rotate(double pivot[], double phi, double theta, double psi) override;
      void translate(double distance[]) override;
   };

   class Region : public TransformableRegion
   {
   public:
      Region(Region* const parent, IMaterialScatterModelT const * const msm, ShapeT const * const shape);
      void removeSubRegion(const RegionBase& subRegion);
      void clearSubRegions();
   };
}

#endif