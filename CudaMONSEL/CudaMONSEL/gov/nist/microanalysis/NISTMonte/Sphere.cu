#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace Sphere
{
   __host__ __device__ Sphere::Sphere(const double center[], double radius) : mRadius(radius)
   {
      memcpy(mCenter, center, sizeof(double) * 3);
   }

   __host__ __device__ bool Sphere::contains(const double pos[]) const
   {
      return Math2::sqr(pos[0] - mCenter[0]) + Math2::sqr(pos[1] - mCenter[1]) + Math2::sqr(pos[2] - mCenter[2]) <= Math2::sqr(mRadius);
   }

   double Sphere::getRadius() const
   {
      return mRadius;
   }

   __host__ __device__ double Sphere::getFirstIntersection(const double pos0[], const double pos1[])
   {
      // Compute the intersection of the line between pos0 and pos1 and the
      // shell of the sphere.
      double d[3];
      Math2::minus3d(pos1, pos0, d);
      double m[3];
      Math2::minus3d(pos0, mCenter, m);
      const double ma2 = -2.0 * Math2::dot3d(d, d);
      const double b = 2.0 * Math2::dot3d(m, d);
      const double c2 = 2.0 * (Math2::dot3d(m, m) - mRadius * mRadius);
      const double f = b * b + ma2 * c2;
      if (f >= 0) {
         double up = (b + ::sqrt(f)) / ma2;
         double un = (b - ::sqrt(f)) / ma2;
         if (up < 0.0)
            up = INFINITY;
         if (un < 0.0)
            un = INFINITY;
         const double res = ::fmin(up, un);
         double s0[3];
         double prd0[3];
         Math2::multiply3d(res, d, prd0);
         Math2::plus3d(m, prd0, s0);
         if (!((res == INFINITY) || (Math2::magnitude3d(s0) - mRadius < ::fmax(1.0e-12, Math2::distance3d(pos0, pos1) * 1.0e-9)))) {
            printf("Sphere::getFirstIntersection wtf: crashed %.10e >= %10e\n", Math2::magnitude3d(s0) - mRadius, ::fmax(1.0e-12, Math2::distance3d(pos0, pos1) * 1.0e-9));
         }
         return res;
      }
      return INFINITY;
   }

   void Sphere::getInitialPoint(int res[]) const
   {
      res[0] = mCenter[0];
      res[1] = mCenter[1];
      res[2] = mCenter[2] - 0.999 * mRadius; // just inside...
   }

   void Sphere::getPointAt(double phi, double theta, double frac, double res[]) const
   {
      res[0] = mCenter[2] + mRadius * frac * ::cos(phi);
      res[1] = mCenter[1] + mRadius * frac * ::sin(phi) * ::sin(theta);
      res[2] = mCenter[0] + mRadius * frac * ::sin(phi) * ::cos(theta);
   }

   // JavaDoc in ITransform
   __host__ __device__ void Sphere::rotate(const double pivot[], double phi, double theta, double psi)
   {
      Transform3D::rotate3d(mCenter, pivot, phi, theta, psi, mCenter);
   }

   // JavaDoc in ITransform
   __host__ __device__ void Sphere::translate(const double distance[])
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

   const double* Sphere::getCenter() const
   {
      return mCenter;
   }

   __host__ __device__ StringT Sphere::toString() const
   {
      return "Sphere[" + amp::to_string(mCenter[0]) + amp::to_string(mCenter[1]) + amp::to_string(mCenter[2]) + ", r=" + amp::to_string(mRadius) + "]";
   }
}