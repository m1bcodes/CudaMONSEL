#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"

#include "gov\nist\nanoscalemetrology\JMONSELutils\NULagrangeInterpolation.cuh"

#include <fstream>

namespace NUTableInterpolation
{
   static std::unordered_map<StringT, const NUTableInterpolationT*> instanceMap;

   /**
   * getInstance - Returns an instance of a RegularTableInterpolation object
   * for the table contained in the named resource.
   *
   * @param tableFileName - A String providing the full path name of the data
   *           file that stores the table to be interpolated.
   */
   const NUTableInterpolation* getInstance(char const * tableFileName)
   {
      printf("NUTableInterpolation* getInstance: %s\n", tableFileName);
      const NUTableInterpolation* uniqueInstance = instanceMap[tableFileName];
      if (uniqueInstance == NULL) {
         uniqueInstance = new NUTableInterpolation(tableFileName); // TODO: fix
         instanceMap[tableFileName] = uniqueInstance;
      }
      return uniqueInstance;
   }

   /**
   * RegularTableInterpolation - Create an interpolation table from the named
   * resource. The table is assumed to be stored in the resource as numbers (in
   * character format) separated by white space. The numbers are in this order:
   * Number of input variables for this table (N), # of values taken by the 1st
   * input variable, the monotonic list of these values, ... (repeated for 2nd,
   * 3rd, up to Nth input variable), then a list of the tabulated values in
   * order with the Nth input variable varying most rapidly, the N-1st next,
   * and so on, with the 1st varying most slowly.
   *
   * @param tableFileName - A String providing the name of the resource (data
   *           file) that stores the table to be interpolated.
   */
   NUTableInterpolation::NUTableInterpolation(char const * tableFileName) : tableFileName(tableFileName), range({ INT_MAX, INT_MIN })
   {
      ReadTable(tableFileName);
   }

   /**
   * Returns the domain (i.e., the interval of valid input values) of the
   * interpolation table. This is an array of double[dim][2] where dim is the
   * number of input variables. domain[i][0] and domain[i][1] are the smallest
   * and largest values of the ith input variable.
   *
   * @return the domain
   */
   MatrixXd NUTableInterpolation::getDomain() const
   {
      // Return a deep copy
      MatrixXd domainCopy(dim, VectorXd(2));
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < 2; j++)
            domainCopy[i][j] = domain[i][j];
      return domainCopy;
   }

   /**
   * Returns a two element array {min,max}, respectively the minimum and
   * maximum output values contained in the table. Note that the interpolated
   * value can fall outside of this range due to overshoot, but usually not by
   * much.
   *
   * @return the range
   */
   VectorXd NUTableInterpolation::getRange() const
   {
      // Return a copy
      VectorXd rangeCopy = {
         range[0],
         range[1]
      };
      return rangeCopy;
   }

   /**
   * interpolate - Interpolates this object's table to determine the value at
   * the supplied input coordinate. If the supplied coordinate lies outside the
   * domain of the table, this method extrapolates. This can very quickly lead
   * to very poor estimates. The calling routine is responsible for checking
   * the input against the domain if extrapolation is to be avoided.
   *
   * @param xval - double[] of length in principle equal to the dimension of
   *           the table. For convenience it is allowed to be greater, in which
   *           case the unnecessary values at the end of the array are ignored.
   * @param order - int The interpolation order, 1 for linear, 3 for cubic,
   *           etc.
   * @return double - The estimated value of the tabulated function at the
   *         supplied coordinate.
   */
   double NUTableInterpolation::interpolate(double xval[], int xvallen, int order) const
   {
      if (xvallen < dim)
         printf("Attempt to interpolate %s at x with %d dimensions", tableFileName.c_str(), dim);

      switch (dim) {
      case 1:
         return NULagrangeInterpolation::d1(table1d.data(), table1d.size(), x[0].data(), x[0].size(), order, xval[0])[0];
      case 2:
         return NULagrangeInterpolation::d2(table2d, x, order, xval, xvallen)[0];
      case 3:
         return NULagrangeInterpolation::d3(table3d, x, order, xval, xvallen)[0];
      case 4:
         return NULagrangeInterpolation::d4(table4d, x, order, xval, xvallen)[0];
      default:
         printf("NUTableInterpolation::interpolate: Table dimensions must be 1<=dim<=4");
      }
   }

   void NUTableInterpolation::ReadTable(char const * tableFileName)
   {
      try {
         std::fstream myfile(tableFileName, std::ios_base::in);
         if (!myfile.good()) throw 0;

         int a;
         while (myfile >> a) {
            dim = (int)a;
            if ((dim < 1) || (dim > 4))
               printf("NUTableInterpolation::ReadTable: Table dimensions must be 1<=dim<=4");
            /*
            * Note: I think I could write a general N-dimension interpolation
            * using techniques similar to Mick Flanagan's PolyCubicSpline
            * algorithm.
            */
            VectorXi nPoints(dim);
            x.resize(dim);

            domain.resize(dim, VectorXd(2));

            for (int i = 0; i < dim; i++) {
               myfile >> a;
               nPoints[i] = a;
               x[i].resize(nPoints[i]);

               for (int j = 0; j < nPoints[i]; j++) {
                  /*
                  * TODO The try/catch below was added to debug a problem that
                  * seems unique to Georg Frase's computer. After the problem is
                  * solved, I can remove it.
                  */
                  myfile >> x[i][j];
               }

               if (x[i][0] < x[i][nPoints[i] - 1]) {
                  domain[i][0] = x[i][0];
                  domain[i][1] = x[i][nPoints[i] - 1];
               }
               else {
                  domain[i][1] = x[i][0];
                  domain[i][0] = x[i][nPoints[i] - 1];
               }
            }

            switch (dim) {
            case 1:
               table1d.resize(nPoints[0]);
               for (int i = 0; i < nPoints[0]; i++) {
                  double tmp;
                  myfile >> tmp;
                  table1d[i] = tmp;
                  if (table1d[i] < range[0])
                     range[0] = table1d[i];
                  else if (table1d[i] > range[1])
                     range[1] = table1d[i];
               }
               break;
            case 2:
               table2d.resize(nPoints[0], VectorXd(nPoints[1]));
               for (int i = 0; i < nPoints[0]; i++)
                  for (int j = 0; j < nPoints[1]; j++) {
                     double tmp;
                     myfile >> tmp;
                     table2d[i][j] = tmp;
                     if (table2d[i][j] < range[0])
                        range[0] = table2d[i][j];
                     else if (table2d[i][j] > range[1])
                        range[1] = table2d[i][j];
                  }
               break;
            case 3:
               table3d.resize(nPoints[0], MatrixXd(nPoints[1], VectorXd(nPoints[2])));
               for (int i = 0; i < nPoints[0]; i++)
                  for (int j = 0; j < nPoints[1]; j++)
                     for (int k = 0; k < nPoints[2]; k++) {
                        double tmp;
                        myfile >> tmp;
                        table3d[i][j][k] = tmp;
                        if (table3d[i][j][k] < range[0])
                           range[0] = table3d[i][j][k];
                        else if (table3d[i][j][k] > range[1])
                           range[1] = table3d[i][j][k];
                     }
               break;
            case 4:
               table4d.resize(nPoints[0], Matrix3DXd(nPoints[1], MatrixXd(nPoints[1], VectorXd(nPoints[2]))));
               for (int i = 0; i < nPoints[0]; i++) {
                  for (int j = 0; j < nPoints[1]; j++) {
                     for (int k = 0; k < nPoints[2]; k++) {
                        for (int m = 0; m < nPoints[3]; m++) {
                           double tmp;
                           myfile >> tmp;
                           table4d[i][j][k][m] = tmp;
                           if (table4d[i][j][k][m] < range[0])
                              range[0] = table4d[i][j][k][m];
                           else if (table4d[i][j][k][m] > range[1])
                              range[1] = table4d[i][j][k][m];
                        }
                     }
                  }
               }
               break;
            }
         }

         myfile.close();
      }
      catch (std::exception&) {
         printf("NUTableInterpolation::ReadTable: failed reading file %s\n", tableFileName);
      }
   }
}
