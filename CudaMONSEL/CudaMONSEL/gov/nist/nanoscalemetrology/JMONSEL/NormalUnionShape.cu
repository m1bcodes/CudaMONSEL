#include "gov\nist\nanoscalemetrology\JMONSEL\NormalUnionShape.cuh"

#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"

namespace NormalUnionShape
{
   __host__ __device__ NormalUnionShape::NormalUnionShape(NormalShapeT& a, NormalShapeT& b) : SumShapeT(&a, &b)
   {  
   }

   __host__ __device__ bool NormalUnionShape::contains(const double pos[]) const
   {
      return SumShapeT::contains(pos);
   }

   __host__ __device__ bool NormalUnionShape::contains(const double pos0[], const double pos1[]) const
   {
      const amp::vector<ShapeT*>& shapes = getShapes();
      for (auto &shape : shapes)
         if (((NormalShapeT*)shape)->contains(pos0, pos1))
            return true;
      return false;
   }

   __host__ __device__ double NormalUnionShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      return getFirstNormal(pos0, pos1)[3];
   }

   __host__ __device__ const double* NormalUnionShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      result[0] = 0.;
      result[1] = 0.;
      result[2] = 0.;
      result[3] = INFINITY;

      const amp::vector<ShapeT*>& shapes = getShapes();
      NormalShape* shapeA = (NormalShape*)shapes.at(0);
      NormalShape* shapeB = (NormalShape*)shapes.at(1);
      int adepth, bdepth;
      double u;

      /*
      * Logic of the algorithm: The algorithm below is a modified form of the
      * "union of dexels" algorithm I designed for a mathematical morphology
      * project. Think of the line parameterized by p = pos0+u*(pos1-pos0) for
      * u from -infinity to infinity. We can characterize the "depth" of parts
      * of this line as 0 if the part is outside both A and B, 1 if it is
      * inside A or B but not both, or 2 if it is inside both. Depth changes
      * only at boundary crossings with shapes A and B, either increasing by 1
      * (for outside-inside) transitions or decreasing by 1 (for inside-outside
      * transitions). (The depth can change by 2 if our line hits A and B
      * boundaries simultaneously.) We are only interested in boundaries for u
      * in (0,1], so we start at u=0 and process boundary crossings in order,
      * updating the depth at each one. If we hit a boundary at which the depth
      * changes from anything to 0 or from 0 to anything, we stop and return
      * this boundary crossing. Depth changes from 1 to 2 or vice versa are
      * "internal boundaries" of the combination, and do not represent
      * boundaries of the union shape, so we keep looking. Alternatively, if we
      * reach u>1 without having found a boundary, we can stop because there
      * are none in (0,1]. The only tricky part is that we don't know when we
      * start at u=0 whether we are inside or outside our object. We could use
      * the contains() methods of the individual shapes, but it is more
      * efficient to wait. Instead we get the first boundary crossing from each
      * shape. (We need this anyway.) If the crossing occurs for u<1, the
      * normal vector is valid, and we can use it to ascertain whether the
      * boundary crossing is inside-out or vice versa. If u>1 the normal vector
      * is not valid, and we have to revert to our fallback plan of using the
      * contains() method.
      */

      // Get 1st A and B intersections and whether we are inside or outside
      const double delta[] = {
         pos1[0] - pos0[0],
         pos1[1] - pos0[1],
         pos1[2] - pos0[2]
      };
      double nva[4];
      memcpy(nva, shapeA->getFirstNormal(pos0, pos1), sizeof(double) * 4);

      if (nva[3] <= 1.)
         adepth = ((delta[0] * nva[0]) + (delta[1] * nva[1]) + (delta[2] * nva[2])) > 0 ? 1 : 0;
      else { // If the crossing is inside-out, then at p0 we are inside.
         // To come back to: What about delta.nva==0?
         if (shapeA->contains(pos0, pos1))
            return result; // If we were inside at p0 and u>1 there
         // can be no inside-out crossing
         adepth = 0;
      }

      double nvb[4];
      memcpy(nvb, shapeB->getFirstNormal(pos0, pos1), sizeof(double) * 4);


      if (nvb[3] <= 1.)
         bdepth = ((delta[0] * nvb[0]) + (delta[1] * nvb[1]) + (delta[2] * nvb[2])) > 0 ? 1 : 0;
      else { // If the crossing is inside-out, then at p0 we are inside.
         // To come back to: What about ?.nvb==0?
         if (shapeB->contains(pos0, pos1))
            return result; // If we were inside at p0 and u>1 there
         // can be no inside-out crossing
         memcpy(result, nva, sizeof(double) * 4);
         return result; // Otherwise there can, and it is at the first A
         // crossing
      }
      int cdepth = adepth + bdepth;

      // See explantion below at first use
      const double EXTRAU = 1.e-10; // This amounts to an Angstrom for a 1
      // meter
      // step.

      for (;;)
         if (nva[3] < nvb[3]) { // shape A provides the first intersection
            if (adepth == cdepth) {
               memcpy(result, nva, sizeof(double) * 4);
               return result; // c toggles from 0 to 1 or vice versa, like A,
               // so this is a boundary
            }

            cdepth = (cdepth == 1 ? 2 : 1); // It wasn't a boundary so we
            // update cdepth

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
            double tmp[] = {
               pos0[0] + (u * delta[0]), // This is pos0+u*delta
               pos0[1] + (u * delta[1]),
               pos0[2] + (u * delta[2])
            };
            memcpy(nva, shapeA->getFirstNormal(tmp, pos1), sizeof(double) * 4); // Find
            // the
            // next
            // one
            // after
            // that

            if (nva[3] < INFINITY)
               nva[3] = (nva[3] * (1. - u)) + u;
            if (nva[3] > 1)
               if (cdepth == bdepth) {
                  memcpy(result, nvb, sizeof(double) * 4);
                  return result;
               }
               else
                  return result;
            adepth = adepth ^ 1; // Toggle depth in A
         }
         else if (nva[3] > nvb[3]) { // Same as above, with A and B roles
            // reversed
            if (bdepth == cdepth) {
               memcpy(result, nvb, sizeof(double) * 4);
               return result; // c toggles from 0 to 1
               // or vice versa, like
               // B, so this is a
               // boundary
            }
            cdepth = (cdepth == 1 ? 2 : 1); // It wasn't a boundary so we
            // update cdepth
            // Get the next intersection in B
            u = nvb[3] + EXTRAU; // Save the distance to our new start
            // point
            double tmp[] = {
               pos0[0] + (u * delta[0]), // This is pos0+u*delta
               pos0[1] + (u * delta[1]),
               pos0[2] + (u * delta[2])
            };
            memcpy(nvb, shapeB->getFirstNormal(tmp, pos1), sizeof(double) * 4);
            // Find
            // the
            // next
            // one
            // after
            // that

            if (nvb[3] < INFINITY)
               nvb[3] = (nvb[3] * (1. - u)) + u;
            if (nvb[3] > 1)
               if (cdepth == adepth) {
                  memcpy(result, nva, sizeof(double) * 4);
                  return result;
               }
               else
                  return result;
            bdepth = bdepth ^ 1; // Toggle depth in B
         }
         else { // Arrive here only in the unlikely event that we
            // simultaneously hit A and B boundaries. Depth changes
            // by 0 or 2
            const int depthchange = (((adepth ^ 1) - adepth) + (bdepth ^ 1)) - bdepth;
            if (depthchange == 0) { // We simultaneously went into one as we
               // went out of the other
               // Update information for both A and B
               // Get the next intersection in both, A first
               u = nva[3] + EXTRAU; // Save the distance to our new
               // start point
               double tmp0[] = {
                  pos0[0] + (u * delta[0]), // This is pos0+u*delta
                  pos0[1] + (u * delta[1]),
                  pos0[2] + (u * delta[2])
               };
               memcpy(nva, shapeB->getFirstNormal(tmp0, pos1), sizeof(double) * 4);

               if (nva[3] < INFINITY)
                  nva[3] = (nva[3] * (1. - u)) + u;

               // Get the next intersection in B
               u = nvb[3] + EXTRAU; // Save the distance to our new
               // start point
               double tmp1[] = {
                  pos0[0] + (u * delta[0]), // This is pos0+u*delta
                  pos0[1] + (u * delta[1]),
                  pos0[2] + (u * delta[2])
               };
               memcpy(nvb, shapeB->getFirstNormal(tmp1, pos1), sizeof(double) * 4); // Find the next one after that
               if (nvb[3] < INFINITY)
                  nvb[3] = (nvb[3] * (1. - u)) + u;

               if (nva[3] > 1)
                  // in A
                  if (cdepth != bdepth) { // remember bdepth changed but
                     // we've not
                     // yet updated it. This really means cdepth ==
                     // bdepth
                     memcpy(result, nvb, sizeof(double) * 4);
                     return result;
                  }
                  else
                     return result;
               if (nvb[3] > 1)
                  // in A
                  if (cdepth != adepth) {// remember adepth changed but
                     // we've not
                     // yet updated it. This really means cdepth ==
                     // adepth
                     memcpy(result, nva, sizeof(double) * 4);
                     return result;
                  }
                  else
                     return result;
               adepth = adepth ^ 1; // Toggle depth in A
               bdepth = bdepth ^ 1; // Toggle depth in B
            }
            else {
               /*
               * Depth went from 0 to 2 or 2 to 0. Either way, this is a
               * boundary. Return average of the two normal vectors. (nva[3]
               * and nvb[3] are the same, so either will do.)
               */
               result[0] = (nva[0] + nvb[0]) / 2.;
               result[0] = (nva[1] + nvb[1]) / 2.,
               result[0] = (nva[2] + nvb[2]) / 2.,
               result[0] = nva[3];
               return result;
            }
         } // End simultaneous boundaries block
   } // End getFirstNormal()

   __host__ __device__ const double* NormalUnionShape::getPreviousNormal() const
   {
      return result;
   }

   __host__ __device__ StringT NormalUnionShape::toString() const
   {
      return "NormalUnionShape";
   }

   __host__ __device__ void NormalUnionShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      SumShapeT::rotate(pivot, phi, theta, psi);
   }

   __host__ __device__ void NormalUnionShape::translate(const double distance[])
   {
      SumShapeT::translate(distance);
   }
}
