#include "gov\nist\nanoscalemetrology\JMONSELutils\NULagrangeInterpolation.cuh"

namespace NULagrangeInterpolation
{
   __host__ __device__ static VectorXi locate(double x, const double xsamp[], int xsamplen, int order)
   {
      /* Use bisection to find the xsamp value nearest x. */

      /* Initialize limits between which the desired index must lie */
      int lowlim = -1;
      int uplim = xsamplen;
      int maxindex = uplim - 1;
      int midpoint;
      int nindex; // nearest index
      bool ascending = xsamp[maxindex] > xsamp[0];

      while ((uplim - lowlim) > 1) {
         midpoint = (uplim + lowlim) >> 1;
         if ((x > xsamp[midpoint]) == ascending)
            lowlim = midpoint;
         else
            uplim = midpoint;
      }

      if (lowlim < 0) {
         int tmp[] = { 0, 0 };
         return VectorXi (tmp, tmp+2);
      }
      else if (uplim > maxindex) {
         int tmp[] = { maxindex - order, maxindex };
         return VectorXi(tmp, tmp + 2);
      }
      if (ascending) {
         if ((x - xsamp[lowlim]) <= (xsamp[uplim] - x))
            nindex = lowlim;
         else
            nindex = uplim;
      }
      else if ((xsamp[lowlim] - x) <= (x - xsamp[uplim]))
         nindex = lowlim;
      else
         nindex = uplim;
      /*
      * At this point lowlim and uplim differ by 1. x is between xsamp[lowlim]
      * and xsamp[uplim]. If order>1 these 2 points are not enough. We need
      * order+1 points for our interpolation. We widen the range by adding
      * points either above or below, adding them in order of proximity to the
      * original interval.
      */
      int firstInd = lowlim; // Initialize
      int lastInd = uplim;

      while (((lastInd - firstInd) < order) && (firstInd > 0) && (lastInd < maxindex))
         if (((((xsamp[lowlim] - xsamp[firstInd - 1]) + xsamp[uplim]) - xsamp[lastInd + 1]) >= 0) == ascending)
            lastInd += 1;
         else
            firstInd -= 1;
      if (lastInd >= maxindex) {
         int tmp[] = { maxindex - order, nindex }; // Extrapolating off high index end of table
         return VectorXi(tmp, tmp + 2);
      }
      else if (firstInd <= 0) {
         int tmp[] = { 0, nindex }; // Extrapolating off the low end
         return VectorXi(tmp, tmp + 2);
      }
      int tmp[] = { firstInd, nindex };
      return VectorXi(tmp, tmp + 2);
   }

   /**
   * neville -- A private utility that performs the actual 1-d interpolation,
   * using Neville's algorithm.
   *
   * @param fsamp - double[] 1D array of function values at the grid points
   * @param xsamp - double[] The values of the independent variable at the grid
   *           points
   * @param location - int[] Location in the table to use, as returned by
   *           locate.
   * @param order - int The order of the interpolation (1 for linear, 3 for
   *           cubic, etc.).
   * @param x - double The x value at which the function is to be estimated.
   * @return double[] - An array of 2 values, the first of which is the
   *         estimate for f[x] and the second is an error estimate.
   */
   __host__ __device__ static VectorXd neville(const double f[], int flen, const double xsamp[], int xsamplen, const int location[], int locationlen, int order, double x)
   {
      int offset = location[0];
      int ns = location[1] - offset; // Nearest grid point

      VectorXd c(f + offset, f + offset + order + 1);
      VectorXd d(f + offset, f + offset + order + 1);

      double y = c[ns--];
      double ho, hp, w, den, dy = 0;
      for (int m = 1; m <= order; m++) {
         for (int i = 0; i <= (order - m); i++) {
            ho = xsamp[i + offset] - x;
            hp = xsamp[i + m + offset] - x;
            w = c[i + 1] - d[i];
            den = ho - hp;
            if (den == 0.)
               printf("NULagrangeInterpolation neville: Identical x values (x = %lf in interpolation table at indices %d and %d.\n", xsamp[i + offset], i + offset, i + offset + m);
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
         }
         dy = ((2 * ns) < (order - 1 - m)) ? c[ns + 1] : d[ns--];
         y += dy;
      }
      VectorXd res(2, 0); res[0] = y; res[1] = dy;
      return res;
   }

   /**
   * d1 - Performs interpolation of the function f(x). The caller supplies a
   * sampling of this function at points xsamp, fsamp. These are arrays (type
   * double[]) of the same length. The values in xsamp must be montonic (not
   * necessarily at equal intervals) without duplicates and the values in fsamp
   * must be the corresponding function values.
   *
   * @param fsamp - double[] 1D array of function values at the grid points
   * @param xsamp - double[] The values of the independent variable at the grid
   *           points
   * @param order - int The order of the interpolation (1 for linear, 3 for
   *           cubic, etc.).
   * @param x - double The x value at which the function value is to be
   *           estimated.
   * @return double[] - An array of 2 values, the first of which is the
   *         estimate for f[x] and the second is an error estimate.
   */
   __host__ __device__ VectorXd d1(const double fsamp[], int fsamplen, const double xsamp[], int xsamplen, int order, double x)
   {
      /*
      * The interpolation is performed by a call to neville, a separate static
      * method.
      */
      if ((order < 1) || (fsamplen < (order + 1))) printf("NULagrangeInterpolation d1: 0 < order <= table.length-1 is required.\n");

      /*
      * index0 is the index of the 1st grid point we'll include in the
      * interpolation. Ideally, we'd like reducedx close to the middle of the
      * set of grid points we choose. If reducedx is close to the beginning or
      * end of the array, this might not be possible.
      */
      VectorXi location = locate(x, xsamp, xsamplen, order);
      return neville(fsamp, fsamplen, xsamp, xsamplen, location.data(), location.size(), order, x);
   }

   /**
   * d2 - Performs interpolation of the 2-d function f(x1,x2). The values of
   * the function on the grid are provided in a double[][] (2-d array). The
   * coordinates of the grid points are assumed to be x1[i] = x0[0]+i*xinc[0]
   * and x2[j] = x0[1]+j*xinc[1] where i is an integer index ranging from 0 to
   * f.length and j is another index ranging from 0 to f[0].length. The grid is
   * specified by providing appropriate x0[] and xinc[] arrays.
   *
   * @param f - double[][] 2-D array of function values at the grid points
   * @param xsamp - double[][] 2 x ? (ragged) array of x values at the grid
   *           points
   * @param order - int The order of the interpolation.
   * @param x - double[] Array of two values, [x1,x2], providing the
   *           coordinates of the point at which the function value is to be
   *           estimated.
   * @return double[] - An array of 2 values, the first of which is the
   *         estimate for f[x] and the second is an error estimate. The error
   *         estimate is based only upon the interpolation in the x1 direction,
   *         however.
   */
   __host__ __device__ VectorXd d2(const MatrixXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen)
   {
      /*
      * N-dimensional interpolation is performed by calling the N-1 dimensional
      * interpolation routine to determine function values at nodes in the Nth
      * dimension, which can then be interpolated via 1-d interpolation.
      */
      if (xlen < 2)
         printf("NULagrangeInterpolation d2: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d2: 0 < order <= table.length-1 is required.\n");

      int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2 */
      VectorXd y(order + 1, 0);
      for (int i = 0; i <= order; i++)
         y[i] = d1(f[index0 + i].data(), f[index0 + i].size(), xsamp[1].data(), xsamp[1].size(), order, x[1])[0];
      /* Make corresponding x array */
      VectorXd xtemp(xsamp[0].data() + index0, xsamp[0].data() + index0 + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }

   /**
   * d3 - Performs interpolation of the 3-d function f(x1,x2,x3). The values of
   * the function on the grid are provided in a double[][][] (3-D array). The
   * coordinates of the grid points are assumed to be x1[i] = x0[0]+i*xinc[0],
   * x2[j] = x0[1]+j*xinc[1], and x3[k] = x0[2]+k*xinc[2] where i, j, and k are
   * integer indices ranging from 0 to 1 less than the lengths of their
   * respective dimensions. The grid is specified by providing appropriate x0[]
   * and xinc[] arrays.
   *
   * @param f - double[][][] 3-D array of function values at the grid points
   * @param xsamp - double[][] 3 x ? (ragged) array of x values at the grid
   *           points
   * @param order - int The order of the interpolation.
   * @param x - double[] Array of 3 values, [x1,x2,x3], providing the
   *           coordinates of the point at which the function value is to be
   *           estimated.
   * @return double[] - An array of 2 values, the first of which is the
   *         estimate for f[x] and the second is an error estimate. The error
   *         estimate is based only upon the interpolation in the x1 direction,
   *         however.
   */
   __host__ __device__ VectorXd d3(const Matrix3DXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen)
   {
      /*
      * N-dimensional interpolation is performed by calling the N-1 dimensional
      * interpolation routine to determine function values at nodes in the Nth
      * dimension, which can then be interpolated via 1-d interpolation.
      */
      if (xlen < 3)
         printf("NULagrangeInterpolation d3: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d3: 0 < order <= table.length-1 is required.\n");

      int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2,x3 */
      VectorXd y(order + 1, 0);
      double reducedx[] = { x[1], x[2] };
      MatrixXd reducedxsamp(2, VectorXd()); reducedxsamp[0] = xsamp[1]; reducedxsamp[1] = xsamp[2];
      for (int i = 0; i <= order; i++)
         y[i] = d2(f[index0 + i], reducedxsamp, order, reducedx, 2)[0];
      /* Make corresponding x array */
      VectorXd xtemp(xsamp[0].data(), xsamp[0].data() + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }

   /**
   * d4 - Performs interpolation of the 4-d function f(x1,x2,x3,x4). The values
   * of the function on the grid are provided in a double[][][][] (4-D array).
   * The coordinates of the grid points are assumed to be x1[i] =
   * x0[0]+i*xinc[0], x2[j] = x0[1]+j*xinc[1], x3[k] = x0[2]+k*xinc[2], x4[m] =
   * x0[3]+m*xinc[3] where i, j, k, and m are integer indices ranging from 0 to
   * 1 less than the lengths of their respective dimensions. The grid is
   * specified by providing appropriate x0[] and xinc[] arrays.
   *
   * @param f - double[][][][] 4-D array of function values at the grid points
   * @param xsamp - double[][] 4 x ? (ragged) array of x values at the grid
   *           points
   * @param order - int The order of the interpolation.
   * @param x - double[] Array of 4 values, [x1,x2,x3,x4], providing the
   *           coordinates of the point at which the function value is to be
   *           estimated.
   * @return double[] - An array of 2 values, the first of which is the
   *         estimate for f[x] and the second is an error estimate. The error
   *         estimate is based only upon the interpolation in the x1 direction,
   *         however.
   */
   __host__ __device__ VectorXd d4(const Matrix4DXd& f, const MatrixXd& xsamp, int order, const double x[], int xlen)
   {
      /*
      * N-dimensional interpolation is performed by calling the N-1 dimensional
      * interpolation routine to determine function values at nodes in the Nth
      * dimension, which can then be interpolated via 1-d interpolation.
      */
      /*
      * I was concerned about the error checks on the following lines and in
      * the similar areas of the lower dimensionality interpolation calls. When
      * d4 repeatedly calls d3, some of the error checking is unnecessary. I
      * mostly use these routines (e.g., in RegularTableInterpolation) in
      * repeated calls with only x differing. On the first call it makes sense
      * to insure the array size is sufficient and interval spacing is nonzero,
      * etc. but on the second call we already know the answer to this.
      * However, by timing some interpolations I observe that removing the
      * error checking completely speeds execution by less than 6%, so the more
      * generally robust error check as currently written is probably worth it.
      */

      if (xlen < 4)
         printf("NULagrangeInterpolation d4: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d4: 0 < order <= table.length-1 is required.\n");

      /*
      * The reduced coordinate is like an index into the array, but it can take
      * on fractional values to indicate positions between the grid points.
      * Obtain reduced 1st coordinate
      */
      int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2,x3,x4 */
      VectorXd y(order + 1, 0);
      double reducedx[] = { x[1], x[2], x[3] };
      MatrixXd reducedxsamp(3, VectorXd()); reducedxsamp[0] = xsamp[1]; reducedxsamp[1] = xsamp[2]; reducedxsamp[2] = xsamp[3];
      for (int i = 0; i <= order; i++)
         y[i] = d3(f[index0 + i], reducedxsamp, order, reducedx, 3)[0];
      /* Make corresponding x array */
      VectorXd xtemp(xsamp[0].data() + index0, xsamp[0].data() + index0 + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }

   __host__ __device__ static VectorXi locate(float x, const float xsamp[], int xsamplen, int order)
   {
      /* Use bisection to find the xsamp value nearest x. */

      /* Initialize limits between which the desired index must lie */
      int lowlim = -1;
      int uplim = xsamplen;
      int maxindex = uplim - 1;
      int midpoint;
      int nindex; // nearest index
      bool ascending = xsamp[maxindex] > xsamp[0];

      while ((uplim - lowlim) > 1) {
         midpoint = (uplim + lowlim) >> 1;
         if ((x > xsamp[midpoint]) == ascending)
            lowlim = midpoint;
         else
            uplim = midpoint;
      }

      if (lowlim < 0) {
         int tmp[] = { 0, 0 };
         return VectorXi(tmp, tmp + 2);
      }
      else if (uplim > maxindex) {
         int tmp[] = { maxindex - order, maxindex };
         return VectorXi(tmp, tmp + 2);
      }
      if (ascending) {
         if ((x - xsamp[lowlim]) <= (xsamp[uplim] - x))
            nindex = lowlim;
         else
            nindex = uplim;
      }
      else if ((xsamp[lowlim] - x) <= (x - xsamp[uplim]))
         nindex = lowlim;
      else
         nindex = uplim;
      /*
      * At this point lowlim and uplim differ by 1. x is between xsamp[lowlim]
      * and xsamp[uplim]. If order>1 these 2 points are not enough. We need
      * order+1 points for our interpolation. We widen the range by adding
      * points either above or below, adding them in order of proximity to the
      * original interval.
      */
      int firstInd = lowlim; // Initialize
      int lastInd = uplim;

      while (((lastInd - firstInd) < order) && (firstInd > 0) && (lastInd < maxindex))
         if (((((xsamp[lowlim] - xsamp[firstInd - 1]) + xsamp[uplim]) - xsamp[lastInd + 1]) >= 0) == ascending)
            lastInd += 1;
         else
            firstInd -= 1;
      if (lastInd >= maxindex) {
         int tmp[] = { maxindex - order, nindex }; // Extrapolating off high index end of table
         return VectorXi(tmp, tmp + 2);
      }
      else if (firstInd <= 0) {
         int tmp[] = { 0, nindex }; // Extrapolating off the low end
         return VectorXi(tmp, tmp + 2);
      }
      int tmp[] = { firstInd, nindex };
      return VectorXi(tmp, tmp + 2);
   }

   __host__ __device__ static VectorXf neville(const float f[], int flen, const float xsamp[], int xsamplen, const int location[], int locationlen, int order, float x)
   {
      const int offset = location[0];
      int ns = location[1] - offset; // Nearest grid point

      VectorXf c(f + offset, f + offset + order + 1);
      VectorXf d(f + offset, f + offset + order + 1);

      double y = c[ns--];
      double ho, hp, w, den, dy = 0;
      for (int m = 1; m <= order; m++) {
         for (int i = 0; i <= (order - m); i++) {
            ho = xsamp[i + offset] - x;
            hp = xsamp[i + m + offset] - x;
            w = c[i + 1] - d[i];
            den = ho - hp;
            if (den == 0.)
               printf("NULagrangeInterpolation neville: Identical x values (x = %lf in interpolation table at indices %d and %d.\n", xsamp[i + offset], i + offset, i + offset + m);
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
         }
         dy = ((2 * ns) < (order - 1 - m)) ? c[ns + 1] : d[ns--];
         y += dy;
      }
      VectorXf res(2, 0); res[0] = y; res[1] = dy;
      return res;
   }

   __host__ __device__ VectorXf d1(const float fsamp[], int fsamplen, const float xsamp[], int xsamplen, int order, float x)
   {
      /*
      * The interpolation is performed by a call to neville, a separate static
      * method.
      */
      if ((order < 1) || (fsamplen < (order + 1))) printf("NULagrangeInterpolation d1: 0 < order <= table.length-1 is required.\n");

      /*
      * index0 is the index of the 1st grid point we'll include in the
      * interpolation. Ideally, we'd like reducedx close to the middle of the
      * set of grid points we choose. If reducedx is close to the beginning or
      * end of the array, this might not be possible.
      */
      VectorXi location = locate(x, xsamp, xsamplen, order);
      return neville(fsamp, fsamplen, xsamp, xsamplen, location.data(), location.size(), order, x);
   }

   __host__ __device__ VectorXf d2(const MatrixXf& f, const MatrixXf& xsamp, int order, const float x[], int xlen)
   {
      /*
      * N-dimensional interpolation is performed by calling the N-1 dimensional
      * interpolation routine to determine function values at nodes in the Nth
      * dimension, which can then be interpolated via 1-d interpolation.
      */
      if (xlen < 2)
         printf("NULagrangeInterpolation d2: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d2: 0 < order <= table.length-1 is required.\n");

      const int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2 */
      VectorXf y(order + 1, 0);
      for (int i = 0; i <= order; i++)
         y[i] = d1(f[index0 + i].data(), f[index0 + i].size(), xsamp[1].data(), xsamp[1].size(), order, x[1])[0];
      /* Make corresponding x array */
      VectorXf xtemp(xsamp[0].data() + index0, xsamp[0].data() + index0 + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }

   __host__ __device__ VectorXf d3(const Matrix3DXf& f, const MatrixXf& xsamp, int order, const float x[], int xlen)
   {
      /*
      * N-dimensional interpolation is performed by calling the N-1 dimensional
      * interpolation routine to determine function values at nodes in the Nth
      * dimension, which can then be interpolated via 1-d interpolation.
      */
      if (xlen < 3)
         printf("NULagrangeInterpolation d3: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d3: 0 < order <= table.length-1 is required.\n");

      const int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2,x3 */
      VectorXf y(order + 1, 0);
      float reducedx[] = { x[1], x[2] };
      MatrixXf reducedxsamp(2, VectorXf()); reducedxsamp[0] = xsamp[1]; reducedxsamp[1] = xsamp[2];
      for (int i = 0; i <= order; i++)
         y[i] = d2(f[index0 + i], reducedxsamp, order, reducedx, 2)[0];
      /* Make corresponding x array */
      VectorXf xtemp(xsamp[0].data(), xsamp[0].data() + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }

   __host__ __device__ VectorXf d4(const Matrix4DXf& f, const MatrixXf& xsamp, int order, const float x[], int xlen)
   {
      if (xlen < 4)
         printf("NULagrangeInterpolation d4: Input array is too short.\n");
      if (f.size() < (order + 1))
         printf("NULagrangeInterpolation d4: 0 < order <= table.length-1 is required.\n");

      const int index0 = locate(x[0], xsamp[0].data(), xsamp[0].size(), order)[0];

      /* Generate and populate an array of function values at x2,x3,x4 */
      VectorXf y(order + 1, 0);
      float reducedx[] = { x[1], x[2], x[3] };
      MatrixXf reducedxsamp(3, VectorXf()); reducedxsamp[0] = xsamp[1]; reducedxsamp[1] = xsamp[2]; reducedxsamp[2] = xsamp[3];
      for (int i = 0; i <= order; i++)
         y[i] = d3(f[index0 + i], reducedxsamp, order, reducedx, 3)[0];
      /* Make corresponding x array */
      VectorXf xtemp(xsamp[0].data() + index0, xsamp[0].data() + index0 + order + 1);

      /* Interpolate these to find the value at x1 */
      return d1(y.data(), y.size(), xtemp.data(), xtemp.size(), order, x[0]);
   }
}
