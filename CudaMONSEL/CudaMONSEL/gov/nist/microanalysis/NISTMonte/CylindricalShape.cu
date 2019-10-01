//package gov.nist.microanalysis.NISTMonte;

#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace CylindricalShape
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static const double EPSILON = 1.0e-40;
#else
   static const double EPSILON = 1.0e-40;
#endif

   __host__ __device__ CylindricalShape::CylindricalShape(const double end0[], const double end1[], double radius)
   {
      memcpy(mEnd0, end0, sizeof(mEnd0[0]) * 3);
      memcpy(mEnd1, end1, sizeof(mEnd1[0]) * 3);
      mRadius2 = radius * radius;
      mDelta[0] = end1[0] - end0[0];
      mDelta[1] = end1[1] - end0[1];
      mDelta[2] = end1[2] - end0[2];

      if (mRadius2 < 1.0e-30) printf("The cylinder radius is unrealistically small.\n");
      mLen2 = Math2::sqr(mDelta[0]) + Math2::sqr(mDelta[1]) + Math2::sqr(mDelta[2]);
      if (mLen2 < 1.0e-30) printf("The cylinder length is unrealistically small.\n");
      mDelta2 = Math2::dot3d(mDelta, mDelta);
   }

   __host__ __device__ CylindricalShape::CylindricalShape(const CylindricalShape& other) :
      mRadius2(other.mRadius2),
      mLen2(other.mLen2),
      mDelta2(other.mDelta2)
   {
      memcpy(mEnd0, other.mEnd0, sizeof(mEnd0[0]) * 3);
      memcpy(mEnd1, other.mEnd1, sizeof(mEnd1[0]) * 3);
      memcpy(mDelta, other.mDelta, sizeof(mDelta[0]) * 3);
   }

   __host__ __device__ double CylindricalShape::closestPointOnAxis(const double p[]) const
   {
      return (mDelta[0] * (p[0] - mEnd0[0]) + mDelta[1] * (p[1] - mEnd0[1]) + mDelta[2] * (p[2] - mEnd0[2])) / mLen2;
   }

   __host__ __device__ double CylindricalShape::distanceSqr(const double p[], double u) const
   {
      return Math2::sqr(p[0] - (mEnd0[0] + u * mDelta[0])) + Math2::sqr(p[1] - (mEnd0[1] + u * mDelta[1])) + Math2::sqr(p[2] - (mEnd0[2] + u * mDelta[2]));
   }

   __host__ __device__ bool CylindricalShape::contains(const double pos[]) const
   {
      // project pos onto the line defined by end0 and end1.
      const double u = closestPointOnAxis(pos);
      // Is this point between end0 and end1 and is pos^2 <= mRadius from the
      // line from end0 to end1?
      return (u >= 0) && (u <= 1.0) && (distanceSqr(pos, u) <= mRadius2);
   }

   __host__ __device__ const double* CylindricalShape::getEnd0() const
   {
      return mEnd0;
   }

   __host__ __device__ const double* CylindricalShape::getEnd1()
   {
      mEnd1[0] = mEnd0[0] + mDelta[0];
      mEnd1[1] = mEnd0[1] + mDelta[1];
      mEnd1[2] = mEnd0[2] + mDelta[2];
      return mEnd1;
   }

   __host__ __device__ static double checkT(double t)
   {
      return t >= 0.0 ? t : INFINITY;
   }

   __host__ __device__ double CylindricalShape::getFirstIntersection(const double sa[], const double sb[])
   {
      if (true) {
         double t0 = INFINITY, t1 = INFINITY, tc = INFINITY;
         double n[3];
         Math2::minus3d(sb, sa, n);
         const double nd = Math2::dot3d(n, mDelta);
         if (nd != 0.0) {
            // Check end cap 0
            double m1[3];
            Math2::minus3d(mEnd0, sa, m1);
            double t = Math2::dot3d(mDelta, m1) / nd;
            if (t > 0.0) {
               double pt[3];
               double ptbtw[3];
               Math2::pointBetween3d(sa, sb, t, ptbtw);
               Math2::minus3d(ptbtw, mEnd0, pt);
               if (Math2::dot3d(pt, pt) < mRadius2)
                  t0 = t;
            }
            // Check end cap 1
            double end1[3];
            Math2::plus3d(mEnd0, mDelta, end1);
            double m2[3];
            Math2::minus3d(end1, sa, m2);
            t = Math2::dot3d(mDelta, m2) / nd;
            if (t > 0.0) {
               double pt[3];
               double ptbtw[3];
               Math2::pointBetween3d(sa, sb, t, ptbtw);
               Math2::minus3d(ptbtw, end1, pt);
               if (Math2::dot3d(pt, pt) < mRadius2)
                  t1 = t;
            }
         }
         const double a = mDelta2 * Math2::dot3d(n, n) - nd * nd;
         if (::fabs(a) > EPSILON) {
            double m[3];
            Math2::minus3d(sa, mEnd0, m);
            const double mn = Math2::dot3d(m, n);
            const double b = mDelta2 * mn - nd * Math2::dot3d(m, mDelta);
            const double md = Math2::dot3d(m, mDelta);
            // Consider the side of the cylinder
            const double c = mDelta2 * (Math2::dot3d(m, m) - mRadius2) - md * md;
            const double discr = b * b - a * c;
            if (discr >= 0.0) {
               const double tm = (-b - ::sqrt(discr)) / a;
               const double tp = (-b + ::sqrt(discr)) / a;
               const double t = ::fmin(tm > 0.0 ? tm : INFINITY, tp > 0.0 ? tp : INFINITY);
               if ((t != INFINITY) && (md + t * nd >= 0.0) && (md + t * nd <= mDelta2))
                  tc = t;
            }
         }
         return ::fmin(t0, ::fmin(t1, tc));
      }
      else {
         double m[3];
         Math2::minus3d(sa, mEnd0, m);
         double n[3];
         Math2::minus3d(sb, sa, n);
         const double md = Math2::dot3d(m, mDelta), nd = Math2::dot3d(n, mDelta), dd = Math2::dot3d(mDelta, mDelta);
         // Segment fully outside end caps...
         if ((md < 0.0) && (md + nd < 0.0))
            return INFINITY;
         if ((md > dd) && (md + nd > dd))
            return INFINITY;
         const double nn = Math2::dot3d(n, n), mn = Math2::dot3d(m, n);
         const double a = dd * nn - nd * nd;
         const double k = Math2::dot3d(m, m) - mRadius2;
         const double c = dd * k - md * md;
         if (::abs(a) < EPSILON) {
            if (md < 0.0)
               return checkT(-mn / nn);
            else if (md > dd)
               return checkT((nd - mn) / nn);
            else
               return 0.0;
         }
         const double b = dd * mn - nd * md;
         const double disc = b * b - a * c;
         if (disc < 0.0)
            return INFINITY;
         double t = (-b - ::sqrt(disc)) / a; // Always a >= 0.0
         if (t < 0.0)
            t = (-b + ::sqrt(disc)) / a;
         double mult0[3];
         Math2::multiply3d(t, n, mult0);
         double p0[3];
         Math2::plus3d(sa, mult0, p0);
         double mult1[3];
         Math2::multiply3d(t, n, mult1);
         double p1[3];
         Math2::plus3d(sa, mult1, p1);

         if (!(::fabs(distanceSqr(p0, closestPointOnAxis(p1) - mRadius2)) < 1.0e-10 * mRadius2)) printf("CylindricalShape::getFirstIntersection: < 1.0e-10 * mRadius2 (%.10e)\n", mRadius2);
         // Check end caps
         if (md + t * nd < 0.0) {
            t = -md / nd;
            return k + 2.0 * t * (mn + t * nn) <= 0.0 ? checkT(t) : INFINITY;
         }
         else if (md + t * nd > dd) {
            t = (dd - md) / nd;
            return k + dd - 2.0 * md + t * (2.0 * (mn - nd) + t * nn) <= 0.0 ? checkT(t) : INFINITY;
         }
         return checkT(t);
      }
   }

   __host__ __device__ void CylindricalShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      Transform3D::rotate3d(mEnd0, pivot, phi, theta, psi, mEnd0);
      Transform3D::rotate3d(mDelta, phi, theta, psi, mDelta);
   }

   __host__ __device__ void CylindricalShape::translate(const double distance[])
   {
      mEnd0[0] += distance[0];
      mEnd0[1] += distance[1];
      mEnd0[2] += distance[2];
   }

   __host__ __device__ double CylindricalShape::getRadius() const
   {
      return ::sqrt(mRadius2);
   }

   double CylindricalShape::getLength() const
   {
      return ::sqrt(mDelta[0] * mDelta[0] + mDelta[1] * mDelta[1] + mDelta[2] * mDelta[2]);
   }

   __host__ __device__ StringT CylindricalShape::toString() const
   {
      StringT res = "CylindricalShape([";
      res += amp::to_string(mEnd0[0]) + "," + amp::to_string(mEnd0[1]) + "," + amp::to_string(mEnd0[2]) + "],[";
      res += amp::to_string(mEnd1[0]) + "," + amp::to_string(mEnd1[1]) + "," + amp::to_string(mEnd1[2]) + "],";
      res += amp::to_string(getRadius()) + ")";
      return res.c_str();
   }

   //public void render(TrajectoryVRML.RenderContext vra, Writer wr)
   //   throws IOException{
   //   final NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
   //   nf.setMaximumFractionDigits(3);
   //   nf.setGroupingUsed(false);
   //   final Color color = vra.getCurrentColor();
   //   final String trStr = nf.format(vra.getTransparency());
   //   final String colorStr = nf.format(color.getRed() / 255.0) + " " + nf.format(color.getGreen() / 255.0) + " "
   //      + nf.format(color.getBlue() / 255.0);
   //   wr.append("\nTransform {\n");
   //   // r is the cross product (1,0,0) x norm(mDelta)
   //   {
   //      final double dm = Math2::magnitude(mDelta);
   //      assert(dm > 0.0);
   //      final double[] r = {
   //         mDelta[2] / dm,
   //         0.0,
   //         -mDelta[0] / dm,
   //      };
   //      final double rm = Math2::magnitude(r);
   //      if (rm > 0.0) { // if rotation required...
   //         nf.setMaximumFractionDigits(5);
   //         double th = ::asin(rm);
   //         if (mDelta[1] < 0.0)
   //            th = ::PI - th;
   //         wr.append(" rotation " + nf.format(r[0] / rm) + " " + nf.format(r[1] / rm) + " " + nf.format(r[2] / rm) + " "
   //            + nf.format(th) + "\n");
   //         nf.setMaximumFractionDigits(3);
   //      }
   //   }
   //   wr.append(" translation " + nf.format((mEnd0[0] + mDelta[0] / 2.0) / TrajectoryVRML.SCALE) + " "
   //      + nf.format((mEnd0[1] + mDelta[1] / 2.0) / TrajectoryVRML.SCALE) + " "
   //      + nf.format((mEnd0[2] + mDelta[2] / 2.0) / TrajectoryVRML.SCALE) + "\n");
   //   wr.append(" children [\n");
   //   wr.append("  Shape {\n");
   //   wr.append("   geometry Cylinder {\n");
   //   wr.append("    radius " + nf.format(getRadius() / TrajectoryVRML.SCALE) + "\n");
   //   wr.append("    height " + nf.format(getLength() / TrajectoryVRML.SCALE) + "\n");
   //   wr.append("    bottom TRUE\n");
   //   wr.append("    side TRUE\n");
   //   wr.append("    top TRUE\n");
   //   wr.append("   }\n");
   //   wr.append("   appearance Appearance {\n");
   //   wr.append("    material Material {\n");
   //   wr.append("     emissiveColor " + colorStr + "\n");
   //   wr.append("     transparency " + trStr + "\n");
   //   wr.append("    }\n");
   //   wr.append("   }\n");
   //   wr.append("  }\n");
   //   wr.append(" ]\n");
   //   wr.append("}");
   //   wr.flush();
   //}
}
