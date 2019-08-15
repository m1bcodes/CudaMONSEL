#ifndef _RREGION_BASE_CUH_
#define _RREGION_BASE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"

#include "Amphibian\Hasher.cuh"

namespace RegionBase
{
   class RegionBase
   {
   public:
      struct PtrHashFcn
      {
         __host__ __device__ inline unsigned int operator() (const RegionBase* ptr) const
         {
            return Hasher::Hash((char *)&ptr, sizeof(ptr));
         }
      };

      struct PtrCompareFcn
      {
         __host__ __device__ inline unsigned int operator() (const RegionBase* lhs, const RegionBase* rhs) const
         {
            return lhs == rhs;
         }
      };

      //typedef std::unordered_set<RegionBase*> RBListT;
      typedef amp::unordered_set<RegionBase*, PtrHashFcn, PtrCompareFcn> RBListT;

   protected:
      TransformableRegion* mParent;
      IMaterialScatterModelT* mScatterModel;
      ShapeT* mShape;
      RBListT mSubRegions;

   public:
      //void updateMaterial(const MaterialT& oldMat, const IMaterialScatterModelT& newMat);
      __host__ __device__ void updateMaterial(const IMaterialScatterModelT& oldMat, IMaterialScatterModelT& newMat);
      __host__ __device__ const MaterialT& getMaterial() const;
      __host__ __device__ IMaterialScatterModelT* getScatterModel() const;
      const RBListT& getSubRegions() const;
      __host__ __device__ void addRegion(RegionBase&); // protected member inaccessble via pointer
      __host__ __device__ const RegionBase* containingSubRegion(const double pos[]) const;
      Element::UnorderedSetT getElements(bool recurse) const;
      __host__ __device__ const RegionBase* findEndOfStep(const double p0[], double p1[]) const; // p1 changes!!!
      __host__ __device__ const ShapeT* getShape() const;
      ShapeT* getShape2();

   //protected:
      double getAtomsPerCubicMeter(const ElementT& el) const;
      __host__ __device__ const RegionBase* getParent() const;
      __host__ __device__ bool isContainingRegion(const RegionBase& searchTarget) const;
      __host__ __device__ StringT toString() const;
   };

   class TransformableRegion : public RegionBase, public ITransformT
   {
   public:
      __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi) override;
      __host__ __device__ void translate(const double distance[]) override;
   };

   class Region : public TransformableRegion
   {
   public:
      __host__ __device__ Region(Region* const parent, IMaterialScatterModelT * const msm, ShapeT * const shape);
      Region(const Region&);
      void removeSubRegion(RegionBase& subRegion);
      void clearSubRegions();
   };
}

#endif