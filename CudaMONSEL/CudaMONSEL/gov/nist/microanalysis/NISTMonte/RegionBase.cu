#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"

namespace RegionBase
{
   static double SMALL_DISP = 1.0e-15; // 1 fm

   RegionBase::RegionBase() : mParent(NULL), mScatterModel(NULL), mShape(NULL)
   {
   }

   RegionBase::RegionBase(const RegionBase& rb) : mParent(rb.mParent), mScatterModel(rb.mScatterModel), mShape(rb.mShape)
   {
      if (&rb == this) return;
   
      //mParent = rb.mParent;
      //mScatterModel = rb.mScatterModel;
      //mShape = rb.mShape;
   
      mSubRegions.clear();
      for (auto sr : rb.mSubRegions) {
         mSubRegions.push_back(sr);
      }
   }

   bool RegionBase::operator==(const RegionBase& rb) const
   {
      return &rb == this;
   }
   
   void RegionBase::updateMaterial(const MaterialT& oldMat, const IMaterialScatterModelT& newMat)
   {
      // Recursively replace all instances of oldMat with newMat
      if (mScatterModel->getMaterial() == oldMat)
         mScatterModel = &newMat;
      for (auto reg : mSubRegions)
         reg->updateMaterial(oldMat, newMat);
   }
   
   void RegionBase::updateMaterial(const IMaterialScatterModelT& oldMat, const IMaterialScatterModelT& newMat)
   {
      // Recursively replace all instances of oldMat with newMat
      if (mScatterModel == &oldMat)
         mScatterModel = &newMat;
      for (auto reg : mSubRegions)
         reg->updateMaterial(oldMat, newMat);
   }
   
   const MaterialT& RegionBase::getMaterial() const
   {
      return mScatterModel->getMaterial();
   }
   
   const IMaterialScatterModelT& RegionBase::getScatterModel() const
   {
      return *mScatterModel;
   }
   
   RegionBase::RBListT RegionBase::getSubRegions() const
   {
      return mSubRegions;
   }

   void RegionBase::addRegion(RegionBase * const rb)
   {
      mSubRegions.push_back(rb);
   }
   
   const RegionBase* RegionBase::containingSubRegion(double pos[]) const
   {
      if (mShape->contains(pos)) {
         for (auto reg : mSubRegions) {
            auto csr = reg->containingSubRegion(pos);
            if (csr != NULL)
               return csr;
         }
         return this;
      }
      return NULL;
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
   
   const RegionBase* RegionBase::findEndOfStep(double p0[], double p1[]) const
   {
      PositionVecT pos0(p0, p0 + 3), pos1(p1, p1 + 3);
      const RegionBase* base = this;
   
      double t = mShape->getFirstIntersection(pos0.data(), pos1.data());
      std::string str(std::string(mShape->toString()) + " " + std::to_string(t));
      if (t < 0.0)  printf("%s\n", str.c_str());
      if ((t <= 1.0) && mParent != NULL)
         base = mParent;
      const RegionBase* res = this;
      for (auto subRegion : mSubRegions) {
         double candidate = subRegion->mShape->getFirstIntersection(pos0.data(), pos1.data());
         if (candidate <= 0.0) printf("%s\n", (std::string(subRegion->mShape->toString()) + " " + std::to_string(candidate)).c_str());
         if ((candidate <= 1.0) && (candidate < t)) {
            t = candidate;
            base = subRegion;
         }
      }
      if (t < 0.0) {
         printf("findEndOfStep\n");
      }
      if (t <= 1.0) {
         PositionVecT delta = Math2::minus(pos1, pos0);
         // Put pos1 exactly on the boundary.
         pos1[0] = pos0[0] + (t * delta[0]);
         pos1[1] = pos0[1] + (t * delta[1]);
         pos1[2] = pos0[2] + (t * delta[2]);
         // Find the region just over the boundary...
         PositionVecT over = Math2::plus(pos1, Math2::multiply(SMALL_DISP, Math2::normalize(delta)));
         while (base != NULL) {
            res = base->containingSubRegion(over.data());
            if (res != NULL)
               return res;
            base = base->mParent;
         }
         return NULL; // new point is nowhere in the chamber
      }
      return res;
   }
   
   const ShapeT& RegionBase::getShape() const
   {
      return *mShape;
   }
   
   double RegionBase::getAtomsPerCubicMeter(const ElementT& el) const
   {
      return getMaterial().atomsPerCubicMeter(el);
   }
   
   const RegionBase& RegionBase::getParent() const
   {
      return *mParent;
   }
   
   bool RegionBase::isContainingRegion(const RegionBase& searchTarget) const
   {
      for (auto sub : mSubRegions)
         if ((sub == &searchTarget) || sub->isContainingRegion(searchTarget))
            return true;
      return false;
   }
   
   char const * RegionBase::toString() const
   {
      std::string str(std::string(mShape->toString()) + " of " + getMaterial().toString());
      return str.c_str();
   }

   void TransformableRegion::rotate(double pivot[], double phi, double theta, double psi)
   {
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
   void TransformableRegion::translate(double distance[])
   {
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

   Region::Region(Region* const parent, IMaterialScatterModelT const * const msm, ShapeT const * const shape)
   {
      mParent = parent;
      mScatterModel = msm;
      mShape = shape;
      if (mParent != NULL) mParent->addRegion(this);
   }

   void Region::removeSubRegion(const RegionBase& subRegion)
   {
      mSubRegions.erase(std::find(mSubRegions.begin(), mSubRegions.end(), &subRegion));
   }

   void Region::clearSubRegions()
   {
      mSubRegions.clear();
   }
}