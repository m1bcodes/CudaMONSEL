#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

namespace RegionBase
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static double SMALL_DISP = 1.0e-15; // 1 fm
#else
   static const double SMALL_DISP = 1.0e-15; // 1 fm
#endif
   
   //void RegionBase::updateMaterial(const MaterialT& oldMat, const IMaterialScatterModelT& newMat)
   //{
   //   // Recursively replace all instances of oldMat with newMat
   //   if (mScatterModel->getMaterial() == oldMat)
   //      mScatterModel = &newMat;
   //   for (auto reg : mSubRegions)
   //      reg->updateMaterial(oldMat, newMat);
   //}
   
   __host__ __device__ void RegionBase::updateMaterial(const IMaterialScatterModelT& oldMat, IMaterialScatterModelT& newMat)
   {
      // Recursively replace all instances of oldMat with newMat
      if (mScatterModel == &oldMat)
         mScatterModel = &newMat;
      for (auto &reg : mSubRegions)
         reg->updateMaterial(oldMat, newMat);
   }
   
   __host__ __device__ const MaterialT& RegionBase::getMaterial() const
   {
      return mScatterModel->getMaterial();
   }
   
   __host__ __device__ IMaterialScatterModelT* RegionBase::getScatterModel() const
   {
      return mScatterModel;
   }
   
   const RegionBase::RBListT& RegionBase::getSubRegions() const
   {
      return mSubRegions;
   }

   __host__ __device__ void RegionBase::addRegion(RegionBase& reg)
   {
      mSubRegions.insert(&reg);
   }
   
   __host__ __device__ const RegionBase* RegionBase::containingSubRegion(const double pos[]) const
   {
      //if (pos[0] != pos[0] || pos[1] != pos[1] || pos[2] != pos[2]) printf("wtf\n");
      if (mShape->contains(pos)) {
         for (auto reg : mSubRegions) {
            const RegionBase* csr = reg->containingSubRegion(pos);
            if (csr) return csr;
         }
         return this;
      }
      return nullptr;
   }
   
   Element::UnorderedSetT RegionBase::getElements(bool recurse) const
   {
      Element::UnorderedSetT res(getMaterial().getElementSet());
      if (recurse) {
         for (auto element : mSubRegions) {
            auto rb = element->getElements(true);
            res.insert(rb.begin(), rb.end());
         }
      }
      return res;
   }
   
   __host__ __device__ const RegionBase* RegionBase::findEndOfStep(const double p0[], double p1[]) const
   {
      if (p1[0] != p1[0] || p1[1] != p1[1] || p1[2] != p1[2])
         printf("RegionBase::findEndOfStep0: (%.5e, %.5e, %.5e)\n", p1[0], p1[1], p1[2]);
      const RegionBase* base = this;
      double t = mShape->getFirstIntersection(p0, p1);
      if (t < 0.0 || t != t) printf("RegionBase::findEndOfStep1: %s %lf\n", mShape->toString().c_str(), t);
      if ((t <= 1.0) && mParent != nullptr)
         base = mParent;

      const RegionBase* res = this;
      for (auto subRegion : mSubRegions) {
         const double candidate = subRegion->mShape->getFirstIntersection(p0, p1);
         if (candidate <= 0.0) printf("RegionBase::findEndOfStep2: %s %lf\n", subRegion->mShape->toString().c_str(), candidate);
         if ((candidate <= 1.0) && (candidate < t)) {
            //printf("%s, %s, %.10e\n", subRegion->toString().c_str(), subRegion->mShape->toString().c_str(), candidate);
            t = candidate;
            base = subRegion;
         }
      }

      if (base == nullptr) printf("RegionBase::findEndOfStep3: base == nullptr 0\n");
      //if (t > -1e-300 && t < 1e-300) t = 0;
      if (t < 0.0) printf("RegionBase::findEndOfStep4: invalid t (%.10e)\n", t);

      if (t <= 1.0) {
         double delta[3];
         Math2::minus3d(p1, p0, delta);
         // Put pos1 exactly on the boundary.
         p1[0] = p0[0] + (t * delta[0]);
         p1[1] = p0[1] + (t * delta[1]);
         p1[2] = p0[2] + (t * delta[2]);
         // Find the region just over the boundary...
         double over[3];
         double prod[3];
         double norm[3];
         Math2::normalize3d(delta, norm);
         Math2::multiply3d(SMALL_DISP, norm, prod);
         Math2::plus3d(p1, prod, over);
         if (over[0] != over[0] || over[1] != over[1] || over[2] != over[2]) printf("RegionBase::findEndOfStep5: wtf\n");
         if (base == nullptr) printf("RegionBase::findEndOfStep6: base == nullptr 1\n");
         while (base != nullptr) {
            res = base->containingSubRegion(over);
            if (res != nullptr)
               return res;
            base = base->mParent;
         }
         return nullptr; // newly created point is nowhere in the chamber
      }
      return res;
   }
   
   __host__ __device__ const ShapeT* RegionBase::getShape() const
   {
      return mShape;
   }

   ShapeT* RegionBase::getShape2()
   {
      return mShape;
   }
   
   double RegionBase::getAtomsPerCubicMeter(const ElementT& el) const
   {
      return getMaterial().atomsPerCubicMeter(el);
   }
   
   __host__ __device__ const RegionBase* RegionBase::getParent() const
   {
      return mParent;
   }
   
   __host__ __device__ bool RegionBase::isContainingRegion(const RegionBase& searchTarget) const
   {
      for (auto sub : mSubRegions)
         if ((sub == &searchTarget) || sub->isContainingRegion(searchTarget))
            return true;
      return false;
   }
   
   __host__ __device__ StringT RegionBase::toString() const
   {
      return mShape->toString() + " of " + getMaterial().toString();
   }

   __host__ __device__ void TransformableRegion::rotate(const double pivot[], double phi, double theta, double psi)
   {
      printf("TransformableRegion::rotate: not supported\n");
      //// check whether we can....
      //if (mShape instanceof ITransform)
      //   for (final Object obj : mSubRegions) {
      //      if (!(obj instanceof ITransform))
      //         throw new EPQFatalException(obj.toString() + " does not support transformation.");
      //   }
      //else
      //   throw new EPQFatalException(mShape.toString() + " does not support transformation.");
      //// then do it...
      //ITransform t = (ITransform)mShape;
      //t.rotate(pivot, phi, theta, psi);
      //for (auto element : mSubRegions) {
      //   t = (ITransform)element;
      //   t.rotate(pivot, phi, theta, psi);
      //}
   }

   // documented in ITransform
   __host__ __device__ void TransformableRegion::translate(const double distance[])
   {
      printf("TransformableRegion::translate: not supported\n");
      //// check whether we can....
      //if (mShape instanceof ITransform)
      //   for (final Object obj : mSubRegions) {
      //      if (!(obj instanceof ITransform))
      //         throw new EPQFatalException(obj.toString() + " does not support transformation.");
      //   }
      //else
      //   throw new EPQFatalException(mShape.toString() + " does not support transformation.");
      //// then do it...
      //ITransform t = (ITransform)mShape;
      //t.translate(distance);
      //for (final Object element : mSubRegions) {
      //   t = (ITransform)element;
      //   t.translate(distance);
      //}
   }

   __host__ __device__ Region::Region(Region* const parent, IMaterialScatterModelT * const msm, ShapeT * const shape)
   {
      mParent = parent;
      mScatterModel = msm;
      mShape = shape;
      if (mParent != nullptr) mParent->addRegion(*this);
      mSubRegions.clear();
   }

   void Region::removeSubRegion(RegionBase& subRegion)
   {
      auto itr = mSubRegions.find(&subRegion);
      if (itr != mSubRegions.end()) mSubRegions.erase(itr);
   }

   void Region::clearSubRegions()
   {
      mSubRegions.clear();
   }
}