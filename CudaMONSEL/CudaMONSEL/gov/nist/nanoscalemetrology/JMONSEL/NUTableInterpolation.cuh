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

      __host__ __device__ const VectorXd& gettable1d() const;
      __host__ __device__ const MatrixXd& gettable2d() const;
      __host__ __device__ const Matrix3DXf& gettable3d() const;
      __host__ __device__ const Matrix4DXd& gettable4d() const;
      __host__ __device__ const MatrixXd& getx() const;
      __host__ __device__ const MatrixXd& getdomain() const;
      __host__ __device__ const VectorXd& getrange() const;
      __host__ __device__ int getdim() const;
      __host__ __device__ StringT gettableFileName() const;

      __device__ void copytable1d(const double*, const unsigned int);
      __device__ void copytable2d(const unsigned int, const double*, const unsigned int);
      __device__ void copytable3d(const unsigned int, const unsigned int, float*, const unsigned int);
      __device__ void copytable4d(const unsigned int, const unsigned int, const unsigned int, const double*, const unsigned int);
      __device__ void copyx(const unsigned int, const double*, const unsigned int);
      __device__ void copydomain(const unsigned int, const double*, const unsigned int);
      __device__ void copyrange(const double*, const unsigned int);
      __device__ void copydim(int);
      __device__ void copytableFileName(const char*);

      __device__ void resizetable2d(const unsigned int);
      __device__ void resizetable3d_0(const unsigned int);
      __device__ void resizetable3d_1(const unsigned int, const unsigned int);
      //__device__ void resizetable3d_2(const unsigned int, const unsigned int, const unsigned int);
      __device__ void resizetable4d_0(const unsigned int);
      __device__ void resizetable4d_1(const unsigned int, const unsigned int);
      __device__ void resizetable4d_2(const unsigned int, const unsigned int, const unsigned int);
      __device__ void resizex(const unsigned int);
      __device__ void resizedomain(const unsigned int);

   private:
      void ReadTable(char const * tableFileName);

      VectorXd table1d;
      MatrixXd table2d;
      //Matrix3DXd table3d;
      Matrix3DXf table3d;
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
      
      __device__ void setInstance(StringT, const NUTableInterpolation*);

   private:
      amp::unordered_map<StringT, const NUTableInterpolationT*, amp::string_cmp, PtrCompareFcn, amp::string_hash, PtrHashFcn> instanceMap;
   };

   extern __host__ __device__ const NUTableInterpolation* getInstance(char const * tableFileName);

   extern __global__ void initFactory();

   extern void copyDataToCuda(char const * tableFileName);
}