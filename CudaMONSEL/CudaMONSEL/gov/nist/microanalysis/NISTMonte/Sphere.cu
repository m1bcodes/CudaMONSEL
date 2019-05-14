#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
//import gov.nist.microanalysis.Utility.Transform3D;

namespace Sphere
{
   Sphere::Sphere(const double center[], double radius) : mRadius(radius), mCenter(center, center+3)
   {
   }

   bool Sphere::contains(const double pos[]) const
   {
      return Math2::sqr(pos[0] - mCenter[0]) + Math2::sqr(pos[1] - mCenter[1]) + Math2::sqr(pos[2] - mCenter[2]) <= Math2::sqr(mRadius);
   }

   double Sphere::getRadius() const
   {
      return mRadius;
   }

   double Sphere::getFirstIntersection(const double pos0[], const double pos1[])
   {
      // Compute the intersection of the line between pos0 and pos1 and the
      // shell of the sphere.
      const PositionVecT d = Math2::minus(PositionVecT(pos1, pos1 + 3), PositionVecT(pos0, pos0 + 3));
      const PositionVecT m = Math2::minus(PositionVecT(pos0, pos0 + 3), mCenter);
      const double ma2 = -2.0 * Math2::dot(d, d);
      const double b = 2.0 * Math2::dot(m, d);
      const double c2 = 2.0 * (Math2::dot(m, m) - mRadius * mRadius);
      const double f = b * b + ma2 * c2;
      if (f >= 0) {
         double up = (b + ::sqrt(f)) / ma2;
         double un = (b - ::sqrt(f)) / ma2;
         if (up < 0.0)
            up = INT_MAX;
         if (un < 0.0)
            un = INT_MAX;
         const double res = ::fmin(up, un);
         if (!((res == INT_MAX) || (Math2::magnitude(Math2::plus(m, Math2::multiply(res, d))) - mRadius < ::fmax(1.0e-12, Math2::distance(PositionVecT(pos0, pos0 + 3), PositionVecT(pos1, pos1 + 3)) * 1.0e-9)))) {
            printf("%s\n", std::to_string(Math2::magnitude(Math2::plus(m, Math2::multiply(res, d))) - mRadius).c_str());
         }
         return res;
      }
      return INT_MAX;
   }

   PositionVecT Sphere::getInitialPoint() const
   {
      double res[3];
      res[0] = mCenter[0];
      res[1] = mCenter[1];
      res[2] = mCenter[2] - 0.999 * mRadius; // just inside...
      return PositionVecT(res, res+3);
   }

   PositionVecT Sphere::getPointAt(double phi, double theta, double frac) const
   {
      PositionVecT res(3, 0);
      res[2] = mCenter[2] + mRadius * frac * ::cos(phi);
      res[1] = mCenter[1] + mRadius * frac * ::sin(phi) * ::sin(theta);
      res[0] = mCenter[0] + mRadius * frac * ::sin(phi) * ::cos(theta);
      return res;
   }

   // JavaDoc in ITransform
   void Sphere::rotate(const double pivot[], double phi, double theta, double psi)
   {
      //mCenter = Transform3D.rotate(mCenter, pivot, phi, theta, psi);
   }

   // JavaDoc in ITransform
   void Sphere::translate(const double distance[])
   {
      mCenter[0] += distance[0];
      mCenter[1] += distance[1];
      mCenter[2] += distance[2];
   }

   //void render(TrajectoryVRML.RenderContext vra, Writer wr)
   //   throws IOException
   //{
   //   final NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
   //   nf.setMaximumFractionDigits(2);
   //   nf.setGroupingUsed(false);
   //   final Color color = vra.getCurrentColor();
   //   final String trStr = nf.format(vra.getTransparency());
   //   final String colorStr = nf.format(color.getRed() / 255.0) + " " + nf.format(color.getGreen() / 255.0) + " "
   //      + nf.format(color.getBlue() / 255.0);
   //   wr.append("Transform {\n");
   //   wr.append(" translation " + nf.format(mCenter[0] / TrajectoryVRML.SCALE) + " "
   //      + nf.format(mCenter[1] / TrajectoryVRML.SCALE) + " " + nf.format(mCenter[2] / TrajectoryVRML.SCALE) + "\n");
   //   wr.append(" children [\n");
   //   wr.append("  Shape {\n");
   //   wr.append("   geometry Sphere { radius " + nf.format(mRadius / TrajectoryVRML.SCALE) + "}\n");
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

   PositionVecT Sphere::getCenter() const
   {
      return mCenter;
   }

   StringT Sphere::toString() const
   {
      return "Sphere[" + std::to_string(mCenter[0]) + std::to_string(mCenter[1]) + std::to_string(mCenter[2]) + ", r=" + std::to_string(mRadius) + "]";
   }
}