#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace Transform3D
{
   __host__ __device__ static MatrixXd rotation(const double phi, const double th, const double psi)
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

   __host__ __device__ static MatrixXd times(const MatrixXd& a, const MatrixXd& b)
   {
      if (a[0].size() != b.size()) {
         printf("matrix multiplication failed: %d, %d\n", a[0].size(), b.size());
         return MatrixXd();
      }
      MatrixXd mult(a.size(), VectorXd(b[0].size(), 0));
      for (int i = 0; i < a.size(); ++i)
         for (int j = 0; j < b[0].size(); ++j)
            for (int k = 0; k < b.size(); ++k)
               mult[i][j] += a[i][k] * b[k][j];
      return mult;
   }

   __host__ __device__ static void times(const double** a, const unsigned int asize0, const unsigned int asize1, const double b[], const unsigned int bsize, double res[])
   {
      if (asize1 != bsize) {
         printf("matrix multiplication failed: %d, %d\n", asize1, bsize);
      }
      memset(res, 0, sizeof(res[0]) * asize0);
      for (int i = 0; i < asize0; ++i) {
         for (int j = 0; j < bsize; ++j)
            res[i] += a[i][j] * b[j];
      }
   }

   __host__ __device__ MatrixXd rot(const double ph, const double th, const double ps)
   {
      return rotation(ph, th, ps);
   }

   // using the x-convention Euler angle
   // http://mathworld.wolfram.com/EulerAngles.html
   __host__ __device__ void rotate3d(const double v[], const double phi, const double th, const double psi, double res[])
   {
      //const MatrixXd& rot = rotation(phi, th, psi);
      double r[3][3];
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
      //MatrixXd m(3, VectorXd(v, v + 3));
      double m[3][3] = {
         { v[0], v[0], v[0] },
         { v[1], v[1], v[1] },
         { v[2], v[2], v[2] },
      };
      //m[0][0] = v[0];
      //m[0][1] = v[0];
      //m[0][2] = v[0];
      //m[1][0] = v[1];
      //m[1][1] = v[1];
      //m[1][2] = v[1];
      //m[2][0] = v[2];
      //m[2][1] = v[2];
      //m[2][2] = v[2];
      //const MatrixXd& prod = times(r, m);
      res[0] = r[0][0] * m[0][0] + r[0][1] * m[1][0] + r[0][2] * m[2][0];
      res[1] = r[1][0] * m[0][0] + r[1][1] * m[1][0] + r[1][2] * m[2][0];
      res[2] = r[2][0] * m[0][0] + r[2][1] * m[1][0] + r[2][2] * m[2][0];
      //res[0] = prod[0][0];
      //res[1] = prod[1][0];
      //res[2] = prod[2][0];

      //printf("%.5e, %.5e, %.5e\n", prod[0][0], prod[0][1], prod[0][2]);
      //printf("%.5e, %.5e, %.5e\n", prod[1][0], prod[1][1], prod[1][2]);
      //printf("%.5e, %.5e, %.5e\n", prod[2][0], prod[2][1], prod[2][2]);

      //res[0] = r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2];
      //res[1] = r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2];
      //res[2] = r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2];

      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n",
      //   rot[0][0] * m[0][0] + rot[0][1] * m[1][0] + rot[0][2] * m[2][0],
      //   rot[1][0] * m[0][0] + rot[1][1] * m[1][0] + rot[1][2] * m[2][0],
      //   rot[2][0] * m[0][0] + rot[2][1] * m[1][0] + rot[2][2] * m[2][0],
      //   rot[0][0] * v[0] + rot[0][1] * v[1] + rot[0][2] * v[2],
      //   rot[1][0] * v[0] + rot[1][1] * v[1] + rot[1][2] * v[2],
      //   rot[2][0] * v[0] + rot[2][1] * v[1] + rot[2][2] * v[2]
      //   );
   }

   __host__ __device__ VectorXd rotate(const double v[], const double phi, const double th, const double psi)
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
      VectorXd ret(3, 0);
      ret[0] = res[0][0];
      ret[1] = res[1][0];
      ret[2] = res[2][0];

      return ret;
   }

   __host__ __device__ VectorXd translate(const double point[], const double distance[], bool negate)
   {
      VectorXd res(3, 0);
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

   __host__ __device__ VectorXd rotate(const double point[], const double pivot[], const double phi, const double theta, const double psi)
   {
      return translate(rotate(translate(point, pivot, true).data(), phi, theta, psi).data(), pivot, false);
   }

   __host__ __device__ void translate3d(const double point[], const double distance[], bool negate, double res[])
   {
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
   }

   __host__ __device__ void rotate3d(const double point[], const double pivot[], const double phi, const double theta, const double psi, double res[])
   {
      double translated1[3];
      translate3d(point, pivot, true, translated1);
      double rotated[3];
      rotate3d(translated1, phi, theta, psi, rotated);
      translate3d(rotated, pivot, false, res);
   }
}