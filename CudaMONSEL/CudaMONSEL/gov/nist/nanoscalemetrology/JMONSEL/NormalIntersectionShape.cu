#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"

namespace NormalIntersectionShape
{
   NormalIntersectionShape::NormalIntersectionShape(NormalShapeT& shapeA, NormalShapeT& shapeB) : shapeA(shapeA), shapeB(shapeB)
   {
   }

   bool NormalIntersectionShape::contains(const double pos[]) const
   {
      return ((ShapeT*)&shapeA)->contains(pos) && ((ShapeT*)&shapeB)->contains(pos);
   }

   bool NormalIntersectionShape::contains(const double pos0[], const double pos1[]) const
   {
      return shapeA.contains(pos0, pos1) && shapeB.contains(pos0, pos1);
   }

   VectorXd NormalIntersectionShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      VectorXd nointersection = { 0., 0., 0., INT_MAX };

      int adepth, bdepth;
      double u;

      /*
      * Logic of the algorithm: This algorithm is adapted from the
      * NormalUnionShape version. See the notes there for a fuller description.
      * The difference here is that for an intersection, we are inside the
      * intersection only when depth = 2. Therefore, instead of depth
      * transitions to and from 0, which were the relevant ones for the union
      * shape, we must here look for transitions to and from depth of 2.
      */

      // Get 1st A and B intersections and whether we are inside or outside
      double delta[] = { pos1[0] - pos0[0], pos1[1] - pos0[1], pos1[2] - pos0[2] };
      auto nva = shapeA.getFirstNormal(pos0, pos1);
      if (nva[3] <= 1.)
         adepth = ((delta[0] * nva[0]) + (delta[1] * nva[1]) + (delta[2] * nva[2])) > 0 ? 1 : 0;
      else { // If the crossing is inside-out, then at p0 we are inside.
         // To come back to: What about delta.nva==0?
         if (!shapeA.contains(pos0, pos1)) // pos0 is outside A, hence
            // outside
            // intersection, and can never enter
            // because there are no A boundary
            // crossings
            return nointersection;
         adepth = 1;
      }

      auto nvb = shapeB.getFirstNormal(pos0, pos1);
      if (nvb[3] <= 1.)
         bdepth = ((delta[0] * nvb[0]) + (delta[1] * nvb[1]) + (delta[2] * nvb[2])) > 0 ? 1 : 0;
      else { // If the crossing is inside-out, then at p0 we are inside.
         // To come back to: What about delta.nva==0?
         if (!shapeB.contains(pos0, pos1)) // pos0 is outside B, hence
            // outside
            // intersection, and can never enter
            // because there are no B boundary
            // crossings
            return nointersection;
         if (adepth == 1) {
            result = nva;
            return nva; // We're inside B. If also inside A then next A
            // crossing is our boundary.
         }
         bdepth = 1;
      }
      int cdepth = adepth + bdepth;

      // See explanation below at first use
      double EXTRAU = 1.e-10; // This amounts to an Angstrom for a 1
      // meter
      // step.

      for (;;)
         if (nva[3] < nvb[3]) { // shape A provides the first intersection
            if (bdepth == 1) { //
               result = nva;
               return nva; // c toggles from 1 to 2
               // or vice versa so this is a boundary
            }
            cdepth = cdepth ^ 1; // bdepth is 0 so cdepth was either 0 or
            // 1.
            // This crossing toggles it.

            // Get the next intersection in A
            /*
            * Round off error can cause problems when we try to calculate a u
            * value for the next intersection after our current one. Our new
            * start position at pos0 + u*delta can, because of roundoff error,
            * be just shy of the boundary. We therefore rediscover the same
            * boundary, at a distance of u = 1.e-16 or so! Unfortunately, since
            * we assume each new boundary crossing toggles inside/outside, this
            * messes us up. To avoid this, we advance the position by a tiny
            * extra amount before looking for the next crossing. This amount
            * must be so small that we don't accidentally skip over any
            * boundaries. Step sizes (pos1-pos0) are likely to be on the order
            * of nanometers, so EXTRAU = 1.e-10 means we're safe as long as our
            * shape boundaries are at least a few times 1.e-19 m apart. Since
            * this is smaller than an atomic nucleus, it seems a safe bet! Even
            * for step sizes up to a meter we're safe for boundaries more than
            * 0.1 nm apart.
            */
            u = nva[3] + EXTRAU; // Save the distance to our new start
            // point
            double tmp[] = { pos0[0] + (u * delta[0]), // This is pos0+u*delta
               pos0[1] + (u * delta[1]),
               pos0[2] + (u * delta[2])
            };
            nva = shapeA.getFirstNormal(tmp, pos1); // Find
            // the
            // next
            // one
            // after
            // that

            if (nva[3] < INT_MAX)
               nva[3] = (nva[3] * (1. - u)) + u;
            if (nva[3] > 1)
               if (cdepth == 0)
                  return nointersection;
               else {
                  result = nvb;
                  return nvb;
               }
               adepth = adepth ^ 1; // Toggle depth in A
         }
         else if (nva[3] > nvb[3]) { // Same as above, with A and B roles
            // reversed
            if (adepth == 1) {//
               result = nvb;
               return nvb; // c toggles from 1 to 2 or vice versa so this
               // is a
               // boundary
            }
            cdepth = cdepth ^ 1; // bdepth is 0 so cdepth was either 0 or
            // 1.
            // This crossing toggles it.

            // Get the next intersection in A
            u = nvb[3] + EXTRAU; // Save the distance to our new start
            // point
            double tmp[] = {
               pos0[0] + (u * delta[0]), // This is pos0+u*delta
               pos0[1] + (u * delta[1]),
               pos0[2] + (u * delta[2])
            };
            nvb = shapeB.getFirstNormal(tmp, pos1); // Find
            // the
            // next
            // one
            // after
            // that

            if (nvb[3] < INT_MAX)
               nvb[3] = (nvb[3] * (1. - u)) + u;
            if (nvb[3] > 1)
               if (cdepth == 0)
                  return nointersection;
               else {
                  result = nva;
                  return nva;
               }
               bdepth = bdepth ^ 1; // Toggle depth in B
         }
         else { // Arrive here only in the unlikely event that we
            // simultaneously hit A and B boundaries. Depth changes
            // by 0 or 2
            int depthchange = (((adepth ^ 1) - adepth) + (bdepth ^ 1)) - bdepth;
            if (depthchange == 0) { // We simultaneously went into one as we
               // went out of the other
               // Update information for both A and B
               // Get the next intersection in both, A first
               u = nva[3] + EXTRAU; // Save the distance to our new
               // start point
               double tmp1[] = {
                  pos0[0] + (u * delta[0]), // This is pos0+u*delta
                  pos0[1] + (u * delta[1]),
                  pos0[2] + (u * delta[2])
               };
               nva = shapeA.getFirstNormal(tmp1, pos1); // Find the next one after that
               if (nva[3] < INT_MAX)
                  nva[3] = (nva[3] * (1. - u)) + u;

               // Get the next intersection in B
               u = nvb[3] + EXTRAU; // Save the distance to our new
               // start point
               double tmp2[] = {
                  pos0[0] + (u * delta[0]), // This is pos0+u*delta
                  pos0[1] + (u * delta[1]),
                  pos0[2] + (u * delta[2])
               };
               nvb = shapeB.getFirstNormal(tmp2, pos1); // Find the next one after that
               if (nvb[3] < INT_MAX)
                  nvb[3] = (nvb[3] * (1. - u)) + u;

               if (nva[3] > 1)
                  // in A
                  if (adepth == 0) {// Remember, we've just had a depth
                     // change in
                     // A and B but have not updated the variables. This
                     // means adepth
                     // is actually 1 and bdepth now actually 0, so next
                     // B
                     // transition is the one we want.
                     result = nvb;
                     return nvb;
                  }
                  else
                     return nointersection;
               if (nvb[3] > 1)
                  // in A
                  if (bdepth == 0) {// Remember, we've just had a depth
                     // change in
                     // A and B but have not updated the variables. This
                     // means bdepth
                     // is actually 1 and adepth now actually 0, so next
                     // A
                     // transition is the one we want.
                     result = nva;
                     return nva;
                  }
                  else
                     return nointersection;
               adepth = adepth ^ 1; // Toggle depth in A
               bdepth = bdepth ^ 1; // Toggle depth in B
            }
            else {
               /*
               * Depth went from 0 to 2 or 2 to 0. Either way, this is a
               * boundary. Return average of the two normal vectors. (nva[3]
               * and nvb[3] are the same, so either will do.)
               */
               result = { (nva[0] + nvb[0]) / 2., (nva[1] + nvb[1]) / 2., (nva[2] + nvb[2]) / 2., nva[3] };
               return result;
            }
         } // End simultaneous boundaries block
   } // End getFirstNormal()

   double NormalIntersectionShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      result = getFirstNormal(pos0, pos1);
      return result[3];
   }

   void NormalIntersectionShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      //if (!(shapeA instanceof ITransform))
      //   throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //((ITransform)shapeA).rotate(pivot, phi, theta, psi);
      //if (!(shapeB instanceof ITransform))
      //   throw new EPQFatalException(shapeB.toString() + " does not support transformation.");
      //((ITransform)shapeB).rotate(pivot, phi, theta, psi);
   }

   void NormalIntersectionShape::translate(const double distance[])
   {
      //if (!(shapeA instanceof ITransform))
      //   throw new EPQFatalException(shapeA.toString() + " does not support transformation.");
      //((ITransform)shapeA).translate(distance);
      //if (!(shapeB instanceof ITransform))
      //   throw new EPQFatalException(shapeB.toString() + " does not support transformation.");
      //((ITransform)shapeB).translate(distance);
   }

   VectorXd NormalIntersectionShape::getPreviousNormal() const
   {
      return result;
   }


   char const * NormalIntersectionShape::toString() const
   {
      return "NormalIntersectionShape";
   }
}