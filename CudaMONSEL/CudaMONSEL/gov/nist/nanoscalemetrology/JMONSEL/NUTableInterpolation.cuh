// file: gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

#include "Amphibian\Hasher.cuh"
#include "Amphibian\unordered_map.cuh"

namespace NUTableInterpolation
{
   struct PtrHashFcn
   {
      __host__ __device__ inline unsigned int operator() (const NUTableInterpolationT* ptr) const
      {
         return Hasher::Hash((char *)&ptr, sizeof(ptr));
      }
   };

   struct PtrCompareFcn
   {
      __host__ __device__ inline unsigned int operator() (const NUTableInterpolationT* lhs, const NUTableInterpolationT* rhs) const
      {
         return lhs == rhs;
      }
   };

   //__host__ __device__ const NUTableInterpolation* getInstance(char const * tableFileName);

   class NUTableInterpolation
   {
   public:
      __host__ __device__ NUTableInterpolation(char const * tableFileName);

      __host__ __device__ double interpolate(double xval[], int xvallen, int order) const;

      __host__ __device__ const MatrixXd& getDomain() const;
      __host__ __device__ const VectorXd& getRange() const;

      __device__ void copytable1d();
      //__device__ void copytable2d();
      //__device__ void copytable3d();
      //__device__ void copytable4d();
      //__device__ void copyx();
      //__device__ void copydomain();
      //__device__ void copyrange();

      __host__ __device__ const VectorXd& gettable1d() const;
      __host__ __device__ const MatrixXd& gettable2d() const;
      __host__ __device__ const Matrix3DXd& gettable3d() const;
      __host__ __device__ const Matrix4DXd& gettable4d() const;
      __host__ __device__ const MatrixXd& getx() const;
      __host__ __device__ const MatrixXd& getdomain() const;
      __host__ __device__ const VectorXd& getrange() const;

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

      //static amp::unordered_map<StringT, const NUTableInterpolationT*, amp::string_cmp, PtrCompareFcn, amp::string_hash, PtrHashFcn> instanceMap;
   };

   class NUTableInterpolationFactory
   {
   public:
      __host__ __device__ NUTableInterpolationFactory();
      __host__ __device__ const NUTableInterpolation* getInstance(char const * tableFileName);

   private:
      amp::unordered_map<StringT, const NUTableInterpolationT*, amp::string_cmp, PtrCompareFcn, amp::string_hash, PtrHashFcn> instanceMap;
   };

   __host__ __device__ const NUTableInterpolation* getInstance(char const * tableFileName);

   __global__ void copytable1d();
   __global__ void initCuda();
}