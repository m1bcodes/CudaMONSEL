#include "gov\nist\nanoscalemetrology\JMONSEL\NormalCylindricalShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace NormalCylindricalShape
{
   __host__ __device__ NormalCylindricalShape::NormalCylindricalShape(const double end0[], const double end1[], double radius) :
      CylindricalShapeT(end0, end1, radius),
      radius2(radius * radius)
   {
      memcpy(this->end0, end0, sizeof(double) * 3);
      Math2::minus3d(end1, end0, axis);
      mLen2 = Math2::dot3d(axis, axis);
      mLen = ::sqrt(mLen2);
      Math2::divide3d(axis, mLen, normalizedaxis);
   }

   __host__ __device__ bool NormalCylindricalShape::contains(const double pos[]) const
   {
      return CylindricalShapeT::contains(pos);
   }

   __host__ __device__ bool NormalCylindricalShape::contains(const double pos0[], const double pos1[]) const
   {
      const double p0c[] = { pos0[0] - end0[0], pos0[1] - end0[1], pos0[2] - end0[2] }; // pos0 - end0
      const double p0cdotn = (p0c[0] * normalizedaxis[0]) + (p0c[1] * normalizedaxis[1]) + (p0c[2] * normalizedaxis[2]);
      const double p0cSquared = (p0c[0] * p0c[0]) + (p0c[1] * p0c[1]) + (p0c[2] * p0c[2]);
      const double r2 = p0cSquared - (p0cdotn * p0cdotn); // Distance of pos0
      // from
      // axis
      const bool on_end0 = (p0cdotn == 0) && (r2 <= radius2);
      const bool on_end1 = (p0cdotn == mLen) && (r2 <= radius2);
      const bool on_cylinder = (r2 == radius2) && (p0cdotn >= 0) && (p0cdotn <= mLen);
      if (on_end0) {
         const double delta[] = {
            pos1[0] - pos0[0],
            pos1[1] - pos0[1],
            pos1[2] - pos0[2]
         };
         // The normal vector on this endcap is -normalizedaxis
         const double deltadotnormalizedaxis = (delta[0] * normalizedaxis[0]) + (delta[1] * normalizedaxis[1]) + (delta[2] * normalizedaxis[2]);
         if (deltadotnormalizedaxis == 0)
            // endcap
            /*
            * If the trajectory lies wholly within the endcap we must resort to
            * an arbitrary assignment.
            */
            if (normalizedaxis[0] == 0) {
               if (normalizedaxis[1] == 0)
                  return normalizedaxis[2] > 0;
               else
                  return normalizedaxis[1] > 0;
            }
            else
               return normalizedaxis[0] > 0;
         return deltadotnormalizedaxis > 0;
      }
      else if (on_end1) {
         const double delta[] = {
            pos1[0] - pos0[0],
            pos1[1] - pos0[1],
            pos1[2] - pos0[2]
         };
         // The normal vector on this endcap is normalizedaxis
         const double deltadotnormalizedaxis = (delta[0] * normalizedaxis[0]) + (delta[1] * normalizedaxis[1]) + (delta[2] * normalizedaxis[2]);
         if (deltadotnormalizedaxis == 0.)
            // endcap
            if (normalizedaxis[0] == 0) {
               if (normalizedaxis[1] == 0)
                  return normalizedaxis[2] > 0;
               else
                  return normalizedaxis[1] > 0;
            }
            else
               return normalizedaxis[0] > 0;
         return deltadotnormalizedaxis < 0;
      }
      else if (on_cylinder) {
         const double delta[] = {
            pos1[0] - pos0[0],
            pos1[1] - pos0[1],
            pos1[2] - pos0[2]
         }; // Vector pointing in direction of
         // trajectory
         const double nvtmp[] = {
            p0c[0] - (p0cdotn * normalizedaxis[0]),
            p0c[1] - (p0cdotn * normalizedaxis[1]),
            p0c[2] - (p0cdotn * normalizedaxis[2])
         }; // An outward
         // pointing vector
         // at pos0
         return ((nvtmp[0] * delta[0]) + (nvtmp[1] * delta[1]) + (nvtmp[2] * delta[2])) < 0;
      }
      // Here if pos0 is not on the boundary. This is the usual case.
      return (r2 < radius2) && (p0cdotn > 0) && (p0cdotn < mLen);
   }

   __host__ __device__ double NormalCylindricalShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      return getFirstNormal(pos0, pos1)[3];
   }

   __host__ __device__ bool NormalCylindricalShape::isNormalShape() const
   {
      return true;
   }

   __host__ __device__ const double* NormalCylindricalShape::getFirstNormal(const double pos0[], const double pos1[])
   {
      nv[0] = 0.;
      nv[1] = 0.;
      nv[2] = 0.;
      nv[3] = INFINITY;

      // Various differences and dot products that we will need:
      /*
      * double[] p0c = Math2.minus(pos0, end0); //pos0 - end0 double p0cSquared
      * = Math2.dot(p0c, p0c); double[] delta = Math2.minus(pos1, pos0); //
      * pos1 - pos0 double p0cdotn = Math2.dot(p0c, normalizedaxis); double
      * deltadotn = Math2.dot(delta, normalizedaxis); double deltadotnSquared =
      * deltadotn * deltadotn; double p0cdotdelta = Math2.dot(p0c, delta);
      * double deltaSquared = Math2.dot(delta, delta); // delta^2
      */

      // Following version avoids function call overhead & shaves about 34%
      // off
      // time!
      const double p0c[] = {
         pos0[0] - end0[0],
         pos0[1] - end0[1],
         pos0[2] - end0[2]
      }; // pos0 - end0
      const double p0cSquared = (p0c[0] * p0c[0]) + (p0c[1] * p0c[1]) + (p0c[2] * p0c[2]);
      const double delta[] = {
         pos1[0] - pos0[0],
         pos1[1] - pos0[1],
         pos1[2] - pos0[2]
      }; // pos1 - pos0
      const double p0cdotn = (p0c[0] * normalizedaxis[0]) + (p0c[1] * normalizedaxis[1]) + (p0c[2] * normalizedaxis[2]);
      const double deltadotn = (delta[0] * normalizedaxis[0]) + (delta[1] * normalizedaxis[1]) + (delta[2] * normalizedaxis[2]);
      const double deltadotnSquared = deltadotn * deltadotn;
      const double p0cdotdelta = (delta[0] * p0c[0]) + (delta[1] * p0c[1]) + (delta[2] * p0c[2]);
      const double deltaSquared = (delta[0] * delta[0]) + (delta[1] * delta[1]) + (delta[2] * delta[2]); // delta^2

      /*
      * There are 4 possible intersections, two through the end caps and two
      * through the cylinder body. We'll number these 1, 2, 3, and 4 to
      * remember which one is the nearest.
      */
      int intersectionnumber = 0;

      double u = 0.;
      double u1 = 0.;
      double u2 = 0.;
      double projection = 0.;// Stores projection of intersection onto
      // cylinder
      // axis
      double savedprojection = 0.;

      // Check end caps
      /*
      * Ignore deltadotaxis=0 as this means trajectory is parallel to the end
      * caps. Such a trajectory either misses the end caps, or if it hits them
      * represents a grazing collision that we can ignore.
      */
      if (deltadotn != 0.) {
         // 1st end cap (end0)
         u1 = -p0cdotn / deltadotn; // Distance to plane containing endcap
         /*
         * rSquared = p0cSquared+2*u1*p0cdotdelta+u1*u1*deltaSquared; is
         * in-plane distance to end0, squared. Condition for a valid endcap
         * intersection is we strike the plane of the endcap within the step
         * and the position of the intersection is within r of the center.
         */
         if ((u1 > 0.) && (u1 <= 1.) && ((p0cSquared + (2 * u1 * p0cdotdelta) + (u1 * u1 * deltaSquared)) <= radius2)) {
            nv[3] = u1;
            intersectionnumber = 1;
         }

         // 2nd end cap (end0+axis)
         u2 = u1 + (mLen / deltadotn); // Distance to plane containing other
         // endcap

         if ((u2 > 0.) && (u2 <= 1.) && (u2 < nv[3])
            && (((p0cSquared + (2 * u2 * p0cdotdelta) + (u2 * u2 * deltaSquared) + mLen2)
            - (2 * mLen * (p0cdotn + (u2 * deltadotn)))) <= radius2)) {
            nv[3] = u2;
            intersectionnumber = 2;
         }
      }

      /*
      * Here we check for intersections in the cylindrical wall. This involves
      * solution of a quadratic equation: soln = (-b/2 +/- sqrt((b/2)^2-a*c))/a
      * In the following lines we construct this solution
      */
      if (!(((u1 < 0) && (u2 < 0)) || ((u1 > 1) && (u2 > 1)))) {
         /*
         * The above "if" statement may or may not improve speed. The logic of
         * the test is this: If u1<0 and u2<0 then pos0 and pos1 are "above"
         * both end caps. If u1>1 && u2>1 they are "below" both end caps. In
         * either case there can be no meaningful intersection with the
         * cylindrical wall either, since any such intersection must be outside
         * the endcaps. Thus, this test may allow us to avoid the following
         * time-consuming calculation, but at the cost of the additional time
         * to conduct the test even in those cases where we must go ahead
         * anyway. Thus, whether the test is a net benefit depends upon the
         * relative number of times we find ourselves in each of these
         * situations. If the cylinder is "encapsulated" in a layer with planes
         * defined by the endcaps (or perhaps even better, inside of a
         * rectangular block), then trajectories outside of this region check
         * for plane intersections (fast) instead of cylindrical ones (slow).
         * They never check the cylinder intersections because the cylinder is
         * not a subregion of the region they occupy. Those inside the
         * subregion always check, as they should, and so spare themselves the
         * cost of this test. With proper sample description the case in which
         * the test is beneficial never arises. Still, the cost of the test
         * appears to be only ~5% in timing trials, so it is probably worth it
         * in case the user doesn't encapsulate his cylinders.
         */

         /*
         * If at least one of the points is outside the cylinder we have to
         * check for intersections. Otherwise we skip it. Probably this test is
         * a net benefit, since for well-designed sample descriptions our
         * trajectory should mainly be called when we are inside the cylinder
         * or at least near it.
         */

         const double r0Squared = p0cSquared - (p0cdotn * p0cdotn); // radius of
         // pos0
         const double r1Squared = (r0Squared + (2. * (p0cdotdelta - (p0cdotn * deltadotn))) + deltaSquared) - deltadotnSquared; // radius
         // of
         // pos1

         if ((r0Squared > radius2) || (r1Squared > radius2)) {
            const double a = deltaSquared - deltadotnSquared;
            const double minusbover2 = (p0cdotn * deltadotn) - p0cdotdelta;
            double term = (minusbover2 * minusbover2) - (a * (r0Squared - radius2));
            if (term >= 0) { // The quadratic has real solutions
               term = ::sqrt(term); // term is now the square root

               u = (minusbover2 + term) / a; // 1st quadratic solution
               if ((u > 0) && (u <= 1) && (u < nv[3])) { // Solution falls
                  // within
                  // step.
                  /*
                  * We check whether the projection of the intersection point
                  * onto the cylinder axis falls between the end caps.
                  */
                  projection = p0cdotn + (u * deltadotn);
                  if ((projection >= 0) && (projection <= mLen)) {
                     nv[3] = u;
                     intersectionnumber = 3;
                     savedprojection = projection; // We'll need this
                     // for normal
                     // vector
                  }
               }

               u = (minusbover2 - term) / a; // 2nd quadratic solution
               if ((u > 0) && (u <= 1) && (u < nv[3])) { // Solution falls
                  // within
                  // step.
                  /*
                  * We check whether the projection of the intersection point
                  * onto the cylinder axis falls between the end caps.
                  */
                  projection = p0cdotn + (u * deltadotn);
                  if ((projection >= 0) && (projection <= mLen)) {
                     nv[3] = u;
                     intersectionnumber = 4;
                     savedprojection = projection; // We'll need this
                     // for normal
                     // vector
                  }
               }
            }
         } // End of "is this a net savings?" if statement
      }

      // Decide which intersection we chose
      switch (intersectionnumber) {
      case 1:
         nv[0] = -normalizedaxis[0];
         nv[1] = -normalizedaxis[1];
         nv[2] = -normalizedaxis[2];
         return nv;
      case 2:
         nv[0] = normalizedaxis[0];
         nv[1] = normalizedaxis[1];
         nv[2] = normalizedaxis[2];
         return nv;
      case 3:
      case 4:
         /*
         * Method 1, by subracting axial component to leave perpendicular
         * one
         */
         const double normalv[] = {
            (p0c[0] + (nv[3] * delta[0])) - (savedprojection * normalizedaxis[0]),
            (p0c[1] + (nv[3] * delta[1])) - (savedprojection * normalizedaxis[1]),
            (p0c[2] + (nv[3] * delta[2])) - (savedprojection * normalizedaxis[2])
         };
         const double normalvmag = ::sqrt((normalv[0] * normalv[0]) + (normalv[1] * normalv[1]) + (normalv[2] * normalv[2]));
         nv[0] = normalv[0] / normalvmag;
         nv[1] = normalv[1] / normalvmag;
         nv[2] = normalv[2] / normalvmag;
         return nv;

         /*
         * Method 2 using cross products Temporarily store p0c+u delta double[]
         * normalv = {p0c[0]+u*delta[0], p0c[1]+u*delta[1],p0c[2]+u*delta[2]};
         * normalv = Math2.normalize(Math2.cross(
         * Math2.cross(normalizedaxis,normalv),normalizedaxis)); result[0] =
         * normalv[0]; result[1] = normalv[1]; result[2] = normalv[2]; return
         * result;
         */
      }
      return nv; // None of the above. There was no intersection.
   }

   __host__ __device__ void NormalCylindricalShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      Transform3D::rotate3d(end0, pivot, phi, theta, psi, end0);
      Transform3D::rotate3d(axis, phi, theta, psi, axis);
      Math2::divide3d(axis, ::sqrt(mLen2), normalizedaxis);
      CylindricalShapeT::rotate(pivot, phi, theta, psi);
   }

   __host__ __device__ void NormalCylindricalShape::translate(const double distance[])
   {
      end0[0] += distance[0];
      end0[1] += distance[1];
      end0[2] += distance[2];
      CylindricalShapeT::translate(distance);
   }

   __host__ __device__ const double* NormalCylindricalShape::getPreviousNormal() const
   {
      return nv;
   }

   __host__ __device__ StringT NormalCylindricalShape::toString() const
   {
      return "NormalCylindricalShape";
   }
}