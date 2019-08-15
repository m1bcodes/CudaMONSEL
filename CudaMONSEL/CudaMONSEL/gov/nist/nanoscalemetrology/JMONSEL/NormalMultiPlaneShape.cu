#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace NormalMultiPlaneShape
{
   __host__ __device__ NormalMultiPlaneShape::NormalMultiPlaneShape()
   {
   }

   __host__ __device__ void NormalMultiPlaneShape::updateCach()
   {
      narray.resize(mPlanes.size(), VectorXd(3));
      carray.resize(mPlanes.size(), VectorXd(3));

      for (int i = 0; i < mPlanes.size(); ++i) {
         narray[i].assign(mPlanes[i]->getNormal(), mPlanes[i]->getNormal() + 3); // This plane's normal vector
         carray[i].assign(mPlanes[i]->getPoint(), mPlanes[i]->getPoint() + 3); // A point in this plane
      }
   }

   __host__ __device__ static bool containsTieBreak(const double normal[])
   {
      if (normal[0] < 0.)
         return false;
      if (normal[0] == 0.) {
         if (normal[1] < 0.)
            return false;
         if (normal[1] == 0.)
            if (normal[2] < 0.)
               return false;
      }
      return true;
   }

   __host__ __device__ const double* NormalMultiPlaneShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      /*
      * Explanation of the algorithm: Consider the line described by
      * pos0+u*(pos1-pos0) for u between -infinity and infinity. Suppose this
      * line intersects the ith plane at ui. Then the part of the line that is
      * inside the half space defined by that plane is an interval, either
      * (-infinity,ui] or [ui,infinity) depending upon whether it is an
      * insidethisplane->outside or outside->insidethisplane transition. The n
      * planes that make up this object define n such intervals. The
      * intersection (in the set theoretic sense) of all of these intervals,
      * call it [umin,umax], is the part of the line that is inside the shape.
      * If umin>0 then our starting position at u=0 is outside of the shape and
      * umin represents the closest boundary crossing (outside to in). If
      * umin<0 then we start inside and umax represents the outward bound
      * crossing. The result can be computed with a single loop through the
      * planes by using two variables, umin and umax, to keep track of the end
      * points of the intersection interval. umax starts at "infinity" and can
      * only decrease. umin starts at -"infinity" and can only increase.
      * Whenever we change an end point we store the index of the current
      * plane. At the end we choose a u value corresponding to the nearest
      * crossing with 0<u<=1 and the corresponding index is used to access the
      * normal vector.
      */

      /*
      * We can optimize by noticing that certain situations signal immediately
      * that no boundary crossing is possible. These situations are: (1) the
      * line is parallel to and outside of any plane. (2) umax < 0 This signals
      * that the entire "inside" is "behind" our line segment, at negative u.
      * (3) umax < umin This signals that we have at least two completely
      * disjoint intervals. There can be no position on this line that is
      * inside all of them. (4) umin > 1 This signals that the entire "inside"
      * is beyond the endpoint of our line segment. If any of these situations
      * arise as we loop through the planes, we immediately abort, returning
      * u>1 and undefined normal vector.
      */

      double umin = -INFINITY; // Starting interval is the
      // whole real line
      double umax = INFINITY;
      double u;
      int minindex = -1; // Stores index of plane responsible for umin
      int maxindex = -1; // Same for umax. Initial values are illegal
      // indices.
      result[0] = 0.;
      result[1] = 0.;
      result[2] = 0.;
      result[3] = INFINITY;
      // Initial value
      // designates
      // no intersection

      const double delta[] = { pos1[0] - pos0[0], pos1[1] - pos0[1], pos1[2] - pos0[2] };
      for (int i = 0; i < mPlanes.size(); i++) {
         double pos0minusc[] = { pos0[0] - carray[i][0], pos0[1] - carray[i][1], pos0[2] - carray[i][2] };
         /*
         * Note significance of the sign of the next two variables numerator<0
         * means pos0 is inside the current plane; numerator>0 means it's
         * outside. denominator<0 means the line segment and the plane normal
         * point in opposite directions. I.e., this intersection is an
         * outside->inside transition. denominator>0 means the opposite.
         * denominator==0 means the trajectory is parallel to the plane.
         */
         const double numerator = (pos0minusc[0] * narray[i][0]) + (pos0minusc[1] * narray[i][1]) + (pos0minusc[2] * narray[i][2]);
         const double denominator = (delta[0] * narray[i][0]) + (delta[1] * narray[i][1]) + (delta[2] * narray[i][2]);
         if (denominator == 0) {
            /*
            * If the trajectory is parallel to the plane there are no
            * intersections. If it starts inside it's always inside. If it
            * starts outside it's always outside. In or Out is determined by
            * the numerator. numerator<0, or =0 with tie break = true, means it
            * is inside. In this case we continue looping, searching for
            * intersections with other planes of this shape. Otherwise, we
            * return u>1.
            */
            if ((numerator < 0) || ((numerator == 0) && containsTieBreak(narray[i].data())))
               continue;
            return result;
         }
         u = -numerator / denominator; // Compute intersection point
         if (denominator > 0) { // This is an insidethisplane->outside
            // transition. It changes umax.
            if (u < umax) {
               if ((u < 0)
                  || (u <= umin)) /*
                                  * If the new umax is < 0 the "inside" is
                                  * behind our line segment If the new umax
                                  * is < umin, this plane's inside and an
                                  * earlier one are disjoint. If umax=umin,
                                  * the trajectory enters and leaves the
                                  * shape at the same point, i.e., it is
                                  * tangent to the surface. Since our shape
                                  * is convex, a line can only be tangent on
                                  * the OUTside, so this counts as a
                                  * non-intersection. In any of these cases,
                                  * abort and return no intersection.
                                  */
                                  return result;
               umax = u;
               maxindex = i; // remember index of this plane
            }
         }
         else if (u > umin) { /*
                              * It changes umin. If the new umin is > 1 the
                              * "inside" is beyond the end of our line
                              * segment. If it is >umax this plane's inside
                              * and an earlier one are disjoint. Return
                              * "no intersection" in either case.
                              */
            if ((u > 1) || (u >= umax))
               return result;
            umin = u;
            minindex = i; // Remember index of this plane
         } // end if
      } // end for

      // When we arrive here [umin,umax] defines the completed intersection
      // interval
      if (umin > 0) { // Our boundary crossing is outside -> inside at umin
         result[3] = umin;
         result[0] = narray[minindex][0];
         result[1] = narray[minindex][1];
         result[2] = narray[minindex][2];
         return result;
      } // Otherwise our starting position was already inside
      if ((umax <= 1) && (umax > 0.)) { // Our boundary crossing is inside ->
         // outside at umax<1
         result[3] = umax;
         result[0] = narray[maxindex][0];
         result[1] = narray[maxindex][1];
         result[2] = narray[maxindex][2];
         return result;
      } // Otherwise the entire pos0, pos1 interval lies inside
      return result; // return "no intersection"
   }

   __host__ __device__ double NormalMultiPlaneShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      return (getFirstNormal(pos0, pos1))[3];
   }

   __host__ __device__ bool NormalMultiPlaneShape::contains(const double pos0[], const double pos1[]) const
   {
      bool didDelta = false;
      double p0cdotn;
      double delta[3];
      // Loop over all planes in the shape
      for (int i = 0; i < mPlanes.size(); i++) {
         const double p0c[] = { pos0[0] - carray[i][0], pos0[1] - carray[i][1], pos0[2] - carray[i][2] };
         p0cdotn = (p0c[0] * narray[i][0]) + (p0c[1] * narray[i][1]) + (p0c[2] * narray[i][2]);
         if (p0cdotn > 0.)
            return false;
         if (p0cdotn == 0.) {
            if (!didDelta) {
               delta[0] = pos1[0] - pos0[0];
               delta[0] = pos1[1] - pos0[1];
               delta[0] = pos1[2] - pos0[2];
               didDelta = true;
            }
            double deltadotn = (delta[0] * narray[i][0]) + (delta[1] * narray[i][1]) + (delta[2] * narray[i][2]);
            if (deltadotn > 0.)
               return false;
            if ((deltadotn == 0.) && !containsTieBreak(narray[i].data()))
               return false;
         }
      }
      return true;
   }

   __host__ __device__ const double* NormalMultiPlaneShape::getPreviousNormal() const
   {
      return result;
   }

   //void addPlane(const double p1[], const double p2[], const double p3[])
   //{
   //   double normal[] = Math2::cross(Math2::minus(p2, p1), Math2::minus(p3, p1));
   //   if (Math2::magnitude(normal) == 0.) printf("addPlane: 3 supplied points must be non-colinear.");
   //   addPlane(normal, p1);
   //}

   //void addPlane(const double normal[], const double point[])
   //{
   //   super.addPlane(normal, point);
   //   updateCach();
   //}

   __host__ __device__ void NormalMultiPlaneShape::addPlane(PlaneT *plane)
   {
      MultiPlaneShapeT::addPlane(plane);
      updateCach();
   }

   __host__ __device__ bool NormalMultiPlaneShape::contains(const double pos[]) const
   {
      if (mPlanes.size() == 1)
         return (mPlanes.at(0))->contains(pos);
      else {
         for (auto &pl : mPlanes)
            if (!pl->contains(pos))
               return false;
         return true;
      }
   }

   int NormalMultiPlaneShape::getNumPlanes() const
   {
      return mPlanes.size();
   }

   //const VectorXd& NormalMultiPlaneShape::getNormal(int index) const
   //{
   //   return narray[index];
   //}

   double NormalMultiPlaneShape::getB(int index) const
   {
      return Math2::dot3d(narray[index].data(), carray[index].data());
   }

   __host__ __device__ void NormalMultiPlaneShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      MultiPlaneShapeT::rotate(pivot, phi, theta, psi);
      updateCach();
   }

   __host__ __device__ void NormalMultiPlaneShape::translate(const double distance[])
   {
      MultiPlaneShapeT::translate(distance);
      updateCach();
   }

   __host__ __device__ StringT NormalMultiPlaneShape::toString() const
   {
      return "NormalMultiPlaneShape";
   }
}
