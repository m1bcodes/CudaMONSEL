//package gov.nist.microanalysis.NISTMonte;

#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace CylindricalShape
{
   static const double EPSILON = 1.0e-40;

   CylindricalShape::CylindricalShape(const double end0[], const double end1[], double radius) : 
      mEnd0(end0, end0 + 3),
      mRadius2(radius * radius)
   {
      mDelta = { end1[0] - end0[0], end1[1] - end0[1], end1[2] - end0[2] };
      if (mRadius2 < 1.0e-30) printf("The cylinder radius is unrealistically small.");
      mLen2 = Math2::sqr(mDelta[0]) + Math2::sqr(mDelta[1]) + Math2::sqr(mDelta[2]);
      if (mLen2 < 1.0e-30) printf("The cylinder length is unrealistically small.");
      mDelta2 = Math2::dot3D(mDelta, mDelta);
   }

   CylindricalShape::CylindricalShape(const CylindricalShape& other) :
      mEnd0(other.mEnd0.begin(), other.mEnd0.end()),
      mDelta(other.mDelta.begin(), other.mDelta.end()),
      mRadius2(other.mRadius2),
      mLen2(other.mLen2),
      mDelta2(other.mDelta2)
   {
   }

   double CylindricalShape::closestPointOnAxis(const double p[]) const
   {
      return (mDelta[0] * (p[0] - mEnd0[0]) + mDelta[1] * (p[1] - mEnd0[1]) + mDelta[2] * (p[2] - mEnd0[2])) / mLen2;
   }

   double CylindricalShape::distanceSqr(const double p[], double u) const
   {
      return Math2::sqr(p[0] - (mEnd0[0] + u * mDelta[0])) + Math2::sqr(p[1] - (mEnd0[1] + u * mDelta[1])) + Math2::sqr(p[2] - (mEnd0[2] + u * mDelta[2]));
   }

   bool CylindricalShape::contains(const double pos[]) const
   {
      // project pos onto the line defined by end0 and end1.
      double u = closestPointOnAxis(pos);
      // Is this point between end0 and end1 and is pos^2 <= mRadius from the
      // line from end0 to end1?
      return (u >= 0) && (u <= 1.0) && (distanceSqr(pos, u) <= mRadius2);
   }

   VectorXd CylindricalShape::getEnd0() const
   {
      return mEnd0;
   }

   VectorXd CylindricalShape::getEnd1() const
   {
      return VectorXd({ mEnd0[0] + mDelta[0], mEnd0[1] + mDelta[1], mEnd0[2] + mDelta[2] });
   }

   static double checkT(double t)
   {
      return t >= 0.0 ? t : INFINITY;
   }

   double CylindricalShape::getFirstIntersection(const double sa[], const double sb[])
   {
      VectorXd saVec(sa, sa + 3), sbVec(sb, sb + 3);
      if (true) {
         double t0 = INFINITY, t1 = INFINITY, tc = INFINITY;
         const VectorXd& n = Math2::minus3d(sb, sa);
         double nd = Math2::dot3D(n, mDelta);
         if (nd != 0.0) {
            // Check end cap 0
            double t = Math2::dot3d(mDelta.data(), Math2::minus3d(mEnd0.data(), sa).data()) / nd;
            if (t > 0.0) {
               const VectorXd& pt = Math2::minus3d(Math2::pointBetween3d(sa, sb, t).data(), mEnd0.data());
               if (Math2::dot3d(pt.data(), pt.data()) < mRadius2)
                  t0 = t;
            }
            // Check end cap 1
            const VectorXd& end1 = Math2::plus(mEnd0, mDelta);
            t = Math2::dot3d(mDelta.data(), Math2::minus3d(end1.data(), sa).data()) / nd;
            if (t > 0.0) {
               const VectorXd& pt = Math2::minus3d(Math2::pointBetween3d(sa, sb, t).data(), end1.data());
               if (Math2::dot3d(pt.data(), pt.data()) < mRadius2)
                  t1 = t;
            }
         }
         double a = mDelta2 * Math2::dot3d(n.data(), n.data()) - nd * nd;
         if (::abs(a) > EPSILON) {
            const VectorXd& m = Math2::minus3d(sa, mEnd0.data());
            double mn = Math2::dot3d(m.data(), n.data());
            double b = mDelta2 * mn - nd * Math2::dot3d(m.data(), mDelta.data());
            double md = Math2::dot3d(m.data(), mDelta.data());
            // Consider the side of the cylinder
            double c = mDelta2 * (Math2::dot3d(m.data(), m.data()) - mRadius2) - md * md;
            double discr = b * b - a * c;
            if (discr >= 0.0) {
               double tm = (-b - ::sqrt(discr)) / a;
               double tp = (-b + ::sqrt(discr)) / a;
               double t = ::fmin(tm > 0.0 ? tm : INFINITY, tp > 0.0 ? tp : INFINITY);
               if ((t != INFINITY) && (md + t * nd >= 0.0) && (md + t * nd <= mDelta2))
                  tc = t;
            }
         }
         return ::fmin(t0, ::fmin(t1, tc));
      }
      else {
         auto m = Math2::minus3d(sa, mEnd0.data()), n = Math2::minus3d(sb, sa);
         double md = Math2::dot3d(m.data(), mDelta.data()), nd = Math2::dot3d(n.data(), mDelta.data()), dd = Math2::dot3d(mDelta.data(), mDelta.data());
         // Segment fully outside end caps...
         if ((md < 0.0) && (md + nd < 0.0))
            return INFINITY;
         if ((md > dd) && (md + nd > dd))
            return INFINITY;
         double nn = Math2::dot3d(n.data(), n.data()), mn = Math2::dot3d(m.data(), n.data());
         double a = dd * nn - nd * nd;
         double k = Math2::dot3d(m.data(), m.data()) - mRadius2;
         double c = dd * k - md * md;
         if (::abs(a) < EPSILON) {
            if (md < 0.0)
               return checkT(-mn / nn);
            else if (md > dd)
               return checkT((nd - mn) / nn);
            else
               return 0.0;
         }
         double b = dd * mn - nd * md;
         double disc = b * b - a * c;
         if (disc < 0.0)
            return INFINITY;
         double t = (-b - ::sqrt(disc)) / a; // Always a >= 0.0
         if (t < 0.0)
            t = (-b + ::sqrt(disc)) / a;
         if (!(::abs(distanceSqr(Math2::plus3d(sa, Math2::multiply3d(t, n.data()).data()).data(), closestPointOnAxis(Math2::plus3d(sa, Math2::multiply3d(t, n.data()).data()).data()) - mRadius2)) < 1.0e-10 * mRadius2)) printf("CylindricalShape::getFirstIntersection: < 1.0e-10 * mRadius2 (%.10e)\n", mRadius2);
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

   void CylindricalShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      mEnd0 = Transform3D::rotate(mEnd0.data(), pivot, phi, theta, psi);
      mDelta = Transform3D::rotate(mDelta.data(), phi, theta, psi);
   }

   void CylindricalShape::translate(const double distance[])
   {
      mEnd0[0] += distance[0];
      mEnd0[1] += distance[1];
      mEnd0[2] += distance[2];
   }

   double CylindricalShape::getRadius() const
   {
      return ::sqrt(mRadius2);
   }

   double CylindricalShape::getLength() const
   {
      return ::sqrt(mDelta[0] * mDelta[0] + mDelta[1] * mDelta[1] + mDelta[2] * mDelta[2]);
   }

   StringT CylindricalShape::toString() const
   {
      StringT res = "Cylinder([";
      double* end1 = getEnd1().data();
      res += std::to_string(getEnd0()[0]) + "," + std::to_string(getEnd0()[1]) + "," + std::to_string(getEnd0()[2]) + "],[";
      res += std::to_string(end1[0]) + "," + std::to_string(end1[1]) + "," + std::to_string(end1[2]) + "],";
      res += std::to_string(getRadius()) + ")";
      return res.c_str();
   }

   ///**
   //* @see gov.nist.microanalysis.NISTMonte.TrajectoryVRML.IRender#render(gov.nist.microanalysis.NISTMonte.TrajectoryVRML.RenderContext,
   //*      java.io.Writer)
   //*/
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
