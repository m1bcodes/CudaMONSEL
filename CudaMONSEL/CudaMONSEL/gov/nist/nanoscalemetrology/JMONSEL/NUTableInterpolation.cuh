// file: gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NUTableInterpolation
{
   class NUTableInterpolation
   {
   public:
      NUTableInterpolation(char const * tableFileName);

      double interpolate(double xval[], int xvallen, int order);

      MatrixXd getDomain() const;
      VectorXd getRange() const;

   private:
      void ReadTable(char const * tableFileName);

      VectorXd table1d;
      MatrixXd table2d;
      Matrix3DXd table3d;
      Matrix4DXd table4d;
      MatrixXd x;
      MatrixXd domain;
      VectorXd range;
      int dim; // dimension of this table
      // int[] nPoints; // Array of length dim with number of points for each x
      // double[] xinc; // Array of length dim with x increment size
      // double[] xmin; // Array of minimum x values
      StringT tableFileName;

      static std::unordered_map<StringT, const NUTableInterpolationT*> instanceMap;
   };
}