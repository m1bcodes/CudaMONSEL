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

      __host__ __device__ float interpolate(float xval[], int xvallen, int order) const;

      __host__ __device__ const VectorXf& gettable1d() const;
      __host__ __device__ const MatrixXf& gettable2d() const;
      __host__ __device__ const Matrix3DXf& gettable3d() const;
      __host__ __device__ const Matrix4DXf& gettable4d() const;
      __host__ __device__ const MatrixXf& getx() const;
      __host__ __device__ const MatrixXf& getdomain() const;
      __host__ __device__ const VectorXf& getrange() const;
      __host__ __device__ int getdim() const;
      __host__ __device__ StringT gettableFileName() const;

      template<typename T>
      __device__ void assigntable1d(T* data, const unsigned int len)
      {
         table1d.set_data(data, len);
      }
      template<typename T>
      __device__ void assigntable2d(const unsigned int i0, T* data, const unsigned int len)
      {
         table2d[i0].set_data(data, len);
      }
      template<typename T>
      __device__ void assigntable3d(const unsigned int i0, const unsigned int i1, T* data, const unsigned int len)
      {
         table3d[i0][i1].set_data(data, len);
      }
      template<typename T>
      __device__ void assigntable4d(const unsigned int i0, const unsigned int i1, const unsigned int d2, T* data, const unsigned int len)
      {
         table4d[i0][i1][d2].set_data(data, len);
      }
      template<typename T>
      __device__ void assignx(const unsigned int i0, T* data, const unsigned int len)
      {
         x[i0].set_data(data, len);
      }
      template<typename T>
      __device__ void assigndomain(const unsigned int i0, T* data, const unsigned int len)
      {
         domain[i0].set_data(data, len);
      }
      template<typename T>
      __device__ void assignrange(T* data, const unsigned int len)
      {
         range.set_data(data, len);
      }
      __device__ void copydim(int);
      __device__ void copytableFileName(const char*);

      __device__ void resizetable2d(const unsigned int);
      __device__ void resizetable3d_0(const unsigned int);
      __device__ void resizetable3d_1(const unsigned int, const unsigned int);
      __device__ void resizetable4d_0(const unsigned int);
      __device__ void resizetable4d_1(const unsigned int, const unsigned int);
      __device__ void resizetable4d_2(const unsigned int, const unsigned int, const unsigned int);
      __device__ void resizex(const unsigned int);
      __device__ void resizedomain(const unsigned int);

   private:
      void ReadTable(char const * tableFileName);

      VectorXf table1d;
      MatrixXf table2d;
      Matrix3DXf table3d;
      Matrix4DXf table4d;
      MatrixXf x;
      MatrixXf domain;
      VectorXf range;
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

   extern void transferDataToCuda(char const * tableFileName);
}