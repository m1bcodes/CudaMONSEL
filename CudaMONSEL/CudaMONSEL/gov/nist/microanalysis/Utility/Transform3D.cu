#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace Transform3D
{
   static MatrixXd rotation(double phi, double th, double psi)
   {
      MatrixXd r(3, VectorXd(3, 0));
      if (phi != 0.0) {
         const double cTh = ::cos(th);
         const double sTh = ::sin(th);
         const double cPhi = ::cos(phi);
         const double sPhi = ::sin(phi);
         const double cPsi = ::cos(psi);
         const double sPsi = ::sin(psi);
         r[0][0] = cPhi * cTh * cPsi - sPhi * sPsi;
         r[0][1] = -sPhi * cTh * cPsi - cPhi * sPsi;
         r[0][2] = sTh * cPsi;
         r[1][0] = sPhi * cPsi + cPhi * cTh * sPsi;
         r[1][1] = -sPhi * cTh * sPsi + cPhi * cPsi;
         r[1][2] = sTh * sPsi;
         r[2][0] = -cPhi * sTh;
         r[2][1] = sTh * sPhi;
         r[2][2] = cTh;
      }
      else {
         // Optimize the common special case phi=0.0
         double cTh = ::cos(th);
         double sTh = ::sin(th);
         double cPsi = ::cos(psi);
         double sPsi = ::sin(psi);
         r[0][0] = cTh * cPsi;
         r[0][1] = -sPsi;
         r[0][2] = sTh * cPsi;
         r[1][0] = cTh * sPsi;
         r[1][1] = cPsi;
         r[1][2] = sTh * sPsi;
         r[2][0] = -sTh;
         r[2][1] = 0.0;
         r[2][2] = cTh;
      }
      return r;
   }

   static MatrixXd times(const MatrixXd& a, const MatrixXd b)
   {
      if (a[0].size() != b.size()) {
         printf("matrix multiplicatino failed: %d, %d\n", a[0].size(), b.size());
         return MatrixXd();
      }
      MatrixXd mult(a.size(), VectorXd(b[0].size()));
      for (int i = 0; i < a.size(); ++i)
         for (int j = 0; j < b[0].size(); ++j)
            for (int k = 0; k < b.size(); ++k)
            {
               mult[i][j] += a[i][k] * b[k][j];
            }
      return mult;
   }

   MatrixXd rot(double ph, double th, double ps)
   {
      return rotation(ph, th, ps);
   }

   /**
   * rotate - Rotate the specified directional vector by the Euler angles phi,
   * th, psi. Rotate the object by phi about the z-axis followed by theta round
   * the y-axis followed by psi around the z-axis.
   *
   * @param vector double[]
   * @param phi double
   * @param th double
   * @param psi double
   * @return double[]
   */
   VectorXd rotate(const double v[], double phi, double th, double psi)
   {
      MatrixXd m(3, VectorXd(v, v + 3));
      m[0][0] = v[0];
      m[0][1] = v[0];
      m[0][2] = v[0];
      m[1][0] = v[1];
      m[1][1] = v[1];
      m[1][2] = v[1];
      m[2][0] = v[2];
      m[2][1] = v[2];
      m[2][2] = v[2];
      MatrixXd res = times(rotation(phi, th, psi), m);
      VectorXd ret = { res[0][0], res[1][0], res[2][0] };
      return ret;
   }

   /**
   * translate - Translate the specified point by the specified distance.
   *
   * @param point double[] - The original point
   * @param distance double[] - The offset
   * @param negate - true to translate -distance or false to translate
   *           +distance
   * @return double[] - A new point translated by the specified amount.
   */
   VectorXd translate(const double point[], const double distance[], bool negate)
   {
      VectorXd res(3);
      if (negate) {
         res[0] = point[0] - distance[0];
         res[1] = point[1] - distance[1];
         res[2] = point[2] - distance[2];
      }
      else {
         res[0] = point[0] + distance[0];
         res[1] = point[1] + distance[1];
         res[2] = point[2] + distance[2];
      }
      return res;
   }

   /**
   * rotate - Rotate the specified point about the center point by the Euler
   * angles specified. Rotate the object around the specified point by phi
   * about the z-axis followed by theta round the y-axis followed by psi around
   * the z-axis.
   *
   * @param point double[] - The point to rotate
   * @param pivot double[] - The pivot point
   * @param phi double - Initial rotation about z
   * @param theta double - Second rotation about y
   * @param psi double - rotation about z
   * @return double[] - The result
   */
   VectorXd rotate(const double point[], const double pivot[], double phi, double theta, double psi)
   {
      return translate(rotate(translate(point, pivot, true).data(), phi, theta, psi).data(), pivot, false);
   }
}