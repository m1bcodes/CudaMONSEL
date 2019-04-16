//#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
//
//void TransformableRegion::rotate(double pivot[], double phi, double theta, double psi)
//{
//   //// check whether we can....
//   //if (mShape instanceof ITransform)
//   //   for (final Object obj : mSubRegions) {
//   //      if (!(obj instanceof ITransform))
//   //         throw new EPQFatalException(obj.toString() + " does not support transformation.");
//   //   }
//   //else
//   //   throw new EPQFatalException(mShape.toString() + " does not support transformation.");
//   //// then do it...
//   //ITransform t = (ITransform)mShape;
//   //t.rotate(pivot, phi, theta, psi);
//   //for (auto element : mSubRegions) {
//   //   t = (ITransform)element;
//   //   t.rotate(pivot, phi, theta, psi);
//   //}
//}
//
//// documented in ITransform
//void TransformableRegion::translate(double distance[])
//{
//   //// check whether we can....
//   //if (mShape instanceof ITransform)
//   //   for (final Object obj : mSubRegions) {
//   //      if (!(obj instanceof ITransform))
//   //         throw new EPQFatalException(obj.toString() + " does not support transformation.");
//   //   }
//   //else
//   //   throw new EPQFatalException(mShape.toString() + " does not support transformation.");
//   //// then do it...
//   //ITransform t = (ITransform)mShape;
//   //t.translate(distance);
//   //for (final Object element : mSubRegions) {
//   //   t = (ITransform)element;
//   //   t.translate(distance);
//   //}
//}
//
//RegionBase::RegionBase(const RegionBase& rb)
//{
//   if (&rb == this) return;
//
//   mParent = rb.mParent;
//   mScatterModel = rb.mScatterModel;
//   mShape = rb.mShape;
//
//   mSubRegions.clear();
//   for (auto sr : rb.mSubRegions) {
//      mSubRegions.push_back(sr);
//   }
//}
//
//bool RegionBase::operator==(const RegionBase& rb)
//{
//   if (&rb == this) {
//      return true;
//   }
//   return false;
//}
//
//void RegionBase::updateMaterial(Material::Material& oldMat, IMaterialScatterModel& newMat)
//{
//   // Recursively replace all instances of oldMat with newMat
//   if (mScatterModel->getMaterial() == oldMat)
//      mScatterModel = &newMat;
//   for (auto reg : mSubRegions)
//      reg->updateMaterial(oldMat, newMat);
//}
//
//void RegionBase::updateMaterial(IMaterialScatterModel& oldMat, IMaterialScatterModel& newMat)
//{
//   // Recursively replace all instances of oldMat with newMat
//   if (mScatterModel == &oldMat)
//      mScatterModel = &newMat;
//   for (auto reg : mSubRegions)
//      reg->updateMaterial(oldMat, newMat);
//}
//
//Material::Material& RegionBase::getMaterial()
//{
//   return mScatterModel->getMaterial();
//}
//
//IMaterialScatterModel& RegionBase::getScatterModel()
//{
//   return *mScatterModel;
//}
//
//RegionBase::RBListT RegionBase::getSubRegions()
//{
//   return mSubRegions;
//}
//
//RegionBase* RegionBase::containingSubRegion(double pos[])
//{
//   if (mShape->contains(pos)) {
//      for (auto reg : mSubRegions) {
//         auto csr = reg->containingSubRegion(pos);
//         if (csr != NULL)
//            return csr;
//      }
//      return this;
//   }
//   return NULL;
//}
//
//Element::UnorderedSetT RegionBase::getElements(bool recurse)
//{
//   Element::UnorderedSetT res(getMaterial().getElementSet());
//   if (recurse) {
//      for (auto element : mSubRegions) {
//         auto rb = element->getElements(true);
//         res.insert(rb.begin(), rb.end());
//      }
//   }
//   return res;
//}
//
//RegionBase* RegionBase::findEndOfStep(double p0[], double p1[])
//{
//   PositionVecT pos0(p0, p0 + 3), pos1(p1, p1 + 3);
//   RegionBase* base = this;
//
//   double t = mShape->getFirstIntersection(pos0.data(), pos1.data());
//   std::string str(std::string(mShape->toString()) + " " + std::to_string(t));
//   if (t < 0.0)  printf("%s\n", str.c_str());
//   if ((t <= 1.0) && mParent != NULL)
//      base = mParent;
//   RegionBase* res = this;
//   for (auto subRegion : mSubRegions) {
//      double candidate = subRegion->mShape->getFirstIntersection(pos0.data(), pos1.data());
//      if (candidate <= 0.0) printf("%s\n", (std::string(subRegion->mShape->toString()) + " " + std::to_string(candidate)).c_str());
//      if ((candidate <= 1.0) && (candidate < t)) {
//         t = candidate;
//         base = subRegion;
//      }
//   }
//   if (t < 0.0) {
//      printf("findEndOfStep\n");
//   }
//   if (t <= 1.0) {
//      PositionVecT delta = Math2::minus(pos1, pos0);
//      // Put pos1 exactly on the boundary.
//      pos1[0] = pos0[0] + (t * delta[0]);
//      pos1[1] = pos0[1] + (t * delta[1]);
//      pos1[2] = pos0[2] + (t * delta[2]);
//      // Find the region just over the boundary...
//      PositionVecT over = Math2::plus(pos1, Math2::multiply(SMALL_DISP, Math2::normalize(delta)));
//      while (base != NULL) {
//         res = base->containingSubRegion(over.data());
//         if (res != NULL)
//            return res;
//         base = base->mParent;
//      }
//      return NULL; // new point is nowhere in the chamber
//   }
//   return res;
//}
//
//Shape& RegionBase::getShape()
//{
//   return *mShape;
//}
//
//double RegionBase::getAtomsPerCubicMeter(Element::Element& el)
//{
//   return getMaterial().atomsPerCubicMeter(el);
//}
//
//RegionBase& RegionBase::getParent()
//{
//   return *mParent;
//}
//
//bool RegionBase::isContainingRegion(RegionBase& searchTarget)
//{
//   for (auto sub : mSubRegions)
//      if ((sub == &searchTarget) || sub->isContainingRegion(searchTarget))
//         return true;
//   return false;
//}
//
//char const * const RegionBase::toString()
//{
//   std::string str(std::string(mShape->toString()) + " of " + getMaterial().toString());
//   return str.c_str();
//}