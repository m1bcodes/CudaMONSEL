#include "gov\nist\nanoscalemetrology\JMONSELutils\ULagrangeInterpolation.cuh"

namespace ULagrangeInterpolation
{
   static VectorXd uNeville(const double f[], int offset, int order, double x)
   {
      int ns = (int)::round(x); // Nearest grid point
      if (ns < 0)
         ns = 0;
      else if (ns > order)
         ns = order;
      double* c = new double[order + 1];
      double* d = new double[order + 1];
      memcpy(c, f + offset, sizeof(double) * (order + 1));
      memcpy(d, f + offset, sizeof(double) * (order + 1));

      double y = c[ns--];
      double ho, hp, w, dy = 0;
      for (int m = 1; m <= order; m++) {
         for (int i = 0; i <= order - m; i++) {
            ho = i - x;
            hp = i + m - x;
            w = c[i + 1] - d[i];
            d[i] = -hp * w / m;
            c[i] = -ho * w / m;
         }
         dy = (2 * ns < (order - 1 - m)) ? c[ns + 1] : d[ns--];
         y += dy;
      }
      delete[] c;
      delete[] d;
      VectorXd res(2, 0); res[0] = y; res[1] = dy;
      return res;
   }

   VectorXd d1(double const * const f, int len, double x0, double xinc, int order, double x)
   {
      if (xinc == 0.)
         printf("ULagrangeInterpolation::d1: Interval spacing must be nonzero.\n");
      if ((order < 1) || (len < order + 1))
         printf("ULagrangeInterpolation::d1: 0 < order <= table.length-1 is required.\n");
      const double reducedx = (x - x0) / xinc;
      int index0 = (int)reducedx - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > len - order - 1)
         index0 = len - order - 1;
      return uNeville(f, index0, order, reducedx - index0);
   }

   VectorXd d1(const VectorXd& f, double x0, double xinc, int order, double x)
   {
      return d1(f.data(), f.size(), x0, xinc, order, x);
   }

   VectorXd d2(double const * const * const f, int flen0, int flen1, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen)
   {
      if ((x0len < 2) || (xinclen < 2) || (xlen < 2))
         printf("ULagrangeInterpolation d2: Input array is too short.\n");
      if (xinc[0] == 0.)
         printf("ULagrangeInterpolation d2: Interval spacing must be nonzero.\n");
      if (flen0 < order + 1)
         printf("ULagrangeInterpolation d2: 0 < order <= table.length-1 is required.\n");
      const double reducedx1 = (x[0] - x0[0]) / xinc[0];
      int index0 = (int)reducedx1 - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > flen0 - order - 1)
         index0 = flen0 - order - 1;
      double* y = new double[order + 1];
      for (int i = 0; i <= order; i++) {
         VectorXd temp = d1(f[index0 + i], flen1, x0[1], xinc[1], order, x[1]);
         y[i] = temp[0];
      }

      VectorXd ret = d1(y, order + 1, x0[0] + index0 * xinc[0], xinc[0], order, x[0]);
      delete[] y;
      return ret;
   }

   VectorXd d2(const MatrixXd& f, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen)
   {
      if ((x0len < 2) || (xinclen < 2) || (xlen < 2))
         printf("ULagrangeInterpolation d2: Input array is too short: %d, %d, %d\n", x0len, xinclen, xlen);
      if (xinc[0] == 0.)
         printf("ULagrangeInterpolation d2: Interval spacing must be nonzero.\n");
      if (f.size() < order + 1)
         printf("ULagrangeInterpolation d2: 0 < order <= table.length-1 is required.\n");
      const double reducedx1 = (x[0] - x0[0]) / xinc[0];
      int index0 = (int)reducedx1 - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > f.size() - order - 1)
         index0 = f.size() - order - 1;
      double* y = new double[order + 1];
      for (int i = 0; i <= order; i++) {
         VectorXd temp = d1(f[index0 + i].data(), f[index0 + i].size(), x0[1], xinc[1], order, x[1]);
         y[i] = temp[0];
      }

      VectorXd ret = d1(y, order + 1, x0[0] + index0 * xinc[0], xinc[0], order, x[0]);
      delete[] y;
      return ret;
   }

   VectorXd d3(double const * const * const * const f, int flen0, int flen1, int flen2, const double x0[], int x0len, const double xinc[], int xinclen, int order, const double x[], int xlen)
   {
      if ((x0len < 3) || (xinclen < 3) || (xlen < 3))
         printf("ULagrangeInterpolation d3: Input array is too short: %d, %d, %d\n", x0len, xinclen, xlen);
      if (xinc[0] == 0.)
         printf("ULagrangeInterpolation d3: Interval spacing must be nonzero.\n");
      if (flen0 < order + 1)
         printf("ULagrangeInterpolation d3: 0 < order <= table.length-1 is required.\n");

      const double reducedx1 = (x[0] - x0[0]) / xinc[0];
      int index0 = (int)reducedx1 - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > flen0 - order - 1)
         index0 = flen0 - order - 1;
      double* y = new double[order + 1];
      const double x0temp[] = { x0[1], x0[2] };
      const double xinctemp[] = { xinc[1], xinc[2] };
      const double xtemp[] = { x[1], x[2] };
      for (int i = 0; i <= order; i++) {
         auto temp = d2(f[index0 + i], flen1, flen2, x0temp, 2, xinctemp, 2, order, xtemp, 2);
         y[i] = temp[0];
      }

      auto ret = d1(y, order + 1, x0[0] + index0 * xinc[0], xinc[0], order, x[0]);
      delete[] y;
      return ret;
   }

   VectorXd d4(double const * const * const * const * const f, int flen0, int flen1, int flen2, int flen3, double x0[], int x0len, double xinc[], int xinclen, int order, double x[], int xlen)
   {
      if ((x0len < 4) || (xinclen < 4) || (xlen < 4))
         printf("ULagrangeInterpolation d4: Input array is too short: %d, %d, %d\n", x0len, xinclen, xlen);
      if (xinc[0] == 0.)
         printf("ULagrangeInterpolation d4: Interval spacing must be nonzero.\n");
      if (flen0 < order + 1)
         printf("ULagrangeInterpolation d4: 0 < order <= table.length-1 is required.\n");

      const double reducedx1 = (x[0] - x0[0]) / xinc[0];
      int index0 = (int)reducedx1 - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > flen0 - order - 1)
         index0 = flen0 - order - 1;
      double* y = new double[order + 1];
      const double x0temp[] = {
         x0[1],
            x0[2],
            x0[3]
      };
      const double xinctemp[] = {
         xinc[1],
            xinc[2],
            xinc[3]
      };
      const double xtemp[] = {
         x[1],
            x[2],
            x[3]
      };
      for (int i = 0; i <= order; i++) {
         auto temp = d3(f[index0 + i], flen1, flen2, flen3, x0temp, 3, xinctemp, 3, order, xtemp, 3);
         y[i] = temp[0];
      }
      auto ret = d1(y, order + 1, x0[0] + index0 * xinc[0], xinc[0], order, x[0]);
      delete[] y;
      return ret;
   }

   __host__ __device__ static VectorXf uNeville(const float f[], const int offset, const int order, const float x)
   {
      int ns = (int)::round(x); // Nearest grid point
      if (ns < 0)
         ns = 0;
      else if (ns > order)
         ns = order;
      float* c = new float[order + 1];
      float* d = new float[order + 1];
      memcpy(c, f + offset, sizeof(c[0]) * (order + 1));
      memcpy(d, f + offset, sizeof(d[0]) * (order + 1));

      float y = c[ns--];
      float ho, hp, w, dy = 0;
      for (int m = 1; m <= order; m++) {
         for (int i = 0; i <= order - m; i++) {
            ho = i - x;
            hp = i + m - x;
            w = c[i + 1] - d[i];
            d[i] = -hp * w / m;
            c[i] = -ho * w / m;
         }
         dy = (2 * ns < (order - 1 - m)) ? c[ns + 1] : d[ns--];
         y += dy;
      }
      delete[] c;
      delete[] d;
      VectorXf res(2, 0); res[0] = y; res[1] = dy;
      return res;
   }

   __host__ __device__ VectorXf d1(float const * const f, const int len, const float x0, const float xinc, const int order, const float x)
   {
      if (xinc == 0.)
         printf("ULagrangeInterpolation::d1: Interval spacing must be nonzero.\n");
      if ((order < 1) || (len < order + 1))
         printf("ULagrangeInterpolation::d1: 0 < order <= table.length-1 is required.\n");
      const float reducedx = (x - x0) / xinc;
      int index0 = (int)reducedx - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > len - order - 1)
         index0 = len - order - 1;
      return uNeville(f, index0, order, reducedx - index0);
   }

   __host__ __device__ VectorXf d1(const VectorXf& f, const float x0, const float xinc, const int order, const float x)
   {
      return d1(f.data(), f.size(), x0, xinc, order, x);
   }

   __host__ __device__ VectorXf d2(const MatrixXf& f, const float x0[], int x0len, const float xinc[], int xinclen, int order, const float x[], int xlen)
   {
      if ((x0len < 2) || (xinclen < 2) || (xlen < 2))
         printf("ULagrangeInterpolation d2: Input array is too short: %d, %d, %d\n", x0len, xinclen, xlen);
      if (xinc[0] == 0.)
         printf("ULagrangeInterpolation d2: Interval spacing must be nonzero.\n");
      if (f.size() < order + 1)
         printf("ULagrangeInterpolation d2: 0 < order <= table.length-1 is required.\n");
      const float reducedx1 = (x[0] - x0[0]) / xinc[0];
      int index0 = (int)reducedx1 - order / 2;
      if (index0 < 0)
         index0 = 0;
      else if (index0 > f.size() - order - 1)
         index0 = f.size() - order - 1;
      float* y = new float[order + 1];
      for (int i = 0; i <= order; i++) {
         y[i] = d1(f[index0 + i].data(), f[index0 + i].size(), x0[1], xinc[1], order, x[1])[0];
      }

      const VectorXf& ret = d1(y, order + 1, x0[0] + index0 * xinc[0], xinc[0], order, x[0]);
      delete[] y;
      return ret;
   }
}
