#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"

#include "gov\nist\nanoscalemetrology\JMONSELutils\NULagrangeInterpolation.cuh"

#include <fstream>

#include "CudaUtil.h"

namespace NUTableInterpolation
{
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
   __host__ __device__ NUTableInterpolation::NUTableInterpolation(char const * tableFileName) :
      table1d(0, 0),
      table2d(0, VectorXf(0, 0)),
      table3d(0, MatrixXf(0, VectorXf(0, 0))),
      table4d(0, Matrix3DXf(0, MatrixXf(0, VectorXf(0, 0)))),
      x(0, VectorXf(0, 0)),
      domain(0, VectorXf(0, 0)),
      range(2, 0),
      tableFileName(tableFileName)
   {
      range[0] = INFINITY;
      range[1] = -INFINITY;

#if (!(defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)))
      ReadTable(tableFileName);
#endif
   }

   __host__ __device__ float NUTableInterpolation::interpolate(float xval[], int xvallen, int order) const
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
         return NAN;
      }
   }

   void NUTableInterpolation::ReadTable(char const * tableFileName)
   {
      printf("void NUTableInterpolation::ReadTable: Reading %s\n", tableFileName);
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
            VectorXi nPoints(dim, 0);
            x.resize(dim);

            domain.resize(dim, VectorXf(2, 0));

            for (int i = 0; i < dim; i++) {
               myfile >> a;
               nPoints[i] = a;
               x[i].resize(nPoints[i]);

               for (int j = 0; j < nPoints[i]; j++) {
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
                  myfile >> table1d[i];
                  if (table1d[i] < range[0])
                     range[0] = table1d[i];
                  else if (table1d[i] > range[1])
                     range[1] = table1d[i];
               }
               break;
            case 2:
               table2d.resize(nPoints[0], VectorXf(nPoints[1], 0));
               for (int i = 0; i < nPoints[0]; i++)
                  for (int j = 0; j < nPoints[1]; j++) {
                     myfile >> table2d[i][j];
                     if (table2d[i][j] < range[0])
                        range[0] = table2d[i][j];
                     else if (table2d[i][j] > range[1])
                        range[1] = table2d[i][j];
                  }
               break;
            case 3:
               table3d.resize(nPoints[0], MatrixXf(nPoints[1], VectorXf(nPoints[2], 0)));
               for (int i = 0; i < nPoints[0]; i++)
                  for (int j = 0; j < nPoints[1]; j++)
                     for (int k = 0; k < nPoints[2]; k++) {
                        myfile >> table3d[i][j][k];
                        if (table3d[i][j][k] < range[0])
                           range[0] = table3d[i][j][k];
                        else if (table3d[i][j][k] > range[1])
                           range[1] = table3d[i][j][k];
                     }
               break;
            case 4:
               table4d.resize(nPoints[0], Matrix3DXf(nPoints[1], MatrixXf(nPoints[1], VectorXf(nPoints[2], 0))));
               for (int i = 0; i < nPoints[0]; i++) {
                  for (int j = 0; j < nPoints[1]; j++) {
                     for (int k = 0; k < nPoints[2]; k++) {
                        for (int m = 0; m < nPoints[3]; m++) {
                           myfile >> table4d[i][j][k][m];
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

   __host__ __device__ const VectorXf& NUTableInterpolation::gettable1d() const
   {
      return table1d;
   }

   __host__ __device__ const MatrixXf& NUTableInterpolation::gettable2d() const
   {
      return table2d;
   }

   __host__ __device__ const Matrix3DXf& NUTableInterpolation::gettable3d() const
   {
      return table3d;
   }

   __host__ __device__ const Matrix4DXf& NUTableInterpolation::gettable4d() const
   {
      return table4d;
   }

   __host__ __device__ const MatrixXf& NUTableInterpolation::getx() const
   {
      return x;
   }

   __host__ __device__ const MatrixXf& NUTableInterpolation::getdomain() const
   {
      return domain;
   }

   __host__ __device__ const VectorXf& NUTableInterpolation::getrange() const
   {
      return range;
   }

   __host__ __device__ int NUTableInterpolation::getdim() const
   {
      return dim;
   }

   __host__ __device__ StringT NUTableInterpolation::gettableFileName() const
   {
      return tableFileName;
   }

   __device__ void NUTableInterpolation::copydim(int d)
   {
      dim = d;
   }

   __device__ void NUTableInterpolation::copytableFileName(const char* data)
   {
      tableFileName = data;
   }

   __device__ void NUTableInterpolation::resizetable2d(const unsigned int d0)
   {
      table2d.resize(d0);
   }

   __device__ void NUTableInterpolation::resizetable3d_0(const unsigned int d0)
   {
      table3d.resize(d0);
   }

   __device__ void NUTableInterpolation::resizetable3d_1(const unsigned int i0, const unsigned int d1)
   {
      table3d[i0].resize(d1);
   }

   __device__ void NUTableInterpolation::resizetable4d_0(const unsigned int d0)
   {
      table4d.resize(d0);
   }

   __device__ void NUTableInterpolation::resizetable4d_1(const unsigned int i0, const unsigned int d1)
   {
      table4d[i0].resize(d1);
   }

   __device__ void NUTableInterpolation::resizetable4d_2(const unsigned int i0, const unsigned int i1, const unsigned int d2)
   {
      table4d[i0][i1].resize(d2);
   }

   __device__ void NUTableInterpolation::resizex(const unsigned int d0)
   {
      x.resize(d0);
   }

   __device__ void NUTableInterpolation::resizedomain(const unsigned int d0)
   {
      domain.resize(d0);
   }

   __host__ __device__ NUTableInterpolationFactory::NUTableInterpolationFactory()
   {
   }

   /**
   * getInstance - Returns an instance of a RegularTableInterpolation object
   * for the table contained in the named resource.
   *
   * @param tableFileName - A String providing the full path name of the data
   *           file that stores the table to be interpolated.
   */
   __host__ __device__ const NUTableInterpolation* NUTableInterpolationFactory::getInstance(char const * tableFileName)
   {
      //printf("NUTableInterpolation* getInstance: %s\n", tableFileName);
      const NUTableInterpolation* uniqueInstance = nullptr;
      StringT key(tableFileName);

      if (!instanceMap.ContainsKey(key)) {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
         printf("NUTableInterpolationFactory::getInstance: need to copy %s on to device\n", tableFileName);
#endif
         uniqueInstance = new NUTableInterpolation(tableFileName);
         instanceMap.Put(key, uniqueInstance);
      }
      uniqueInstance = instanceMap[key];

      if (!uniqueInstance) printf("const NUTableInterpolation* NUTableInterpolationFactory::getInstance: failed\n");
      return uniqueInstance;
   }

   __device__ void NUTableInterpolationFactory::setInstance(StringT k, const NUTableInterpolation* v)
   {
      instanceMap.Put(k, v);
   }

   static NUTableInterpolationFactory Factory;
   __device__ static NUTableInterpolationFactory* d_Factory = nullptr;

   __global__ void initFactory()
   {
      d_Factory = new NUTableInterpolationFactory();
   }

   __host__ __device__ const NUTableInterpolation* getInstance(char const * tableFileName)
   {
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
      return d_Factory->getInstance(tableFileName);
#else
      return Factory.getInstance(tableFileName);
#endif
   }

   __device__ NUTableInterpolation* newOne = nullptr;

   __global__ void initNewNUTableOnDevice(char const * fn)
   {
      newOne = new NUTableInterpolation(fn);
   }

   template<typename T>
   __global__ void assigntable1d(T* data, const unsigned int len)
   {
      newOne->assigntable1d(data, len);
   }

   template<typename T>
   __global__ void assigntable2d(const unsigned int i0, T* data, const unsigned int len)
   {
      newOne->assigntable2d(i0, data, len);
   }

   template<typename T>
   __global__ void assigntable3d(const unsigned int i0, const unsigned int i1, T* data, const unsigned int len)
   {
      newOne->assigntable3d(i0, i1, data, len);
   }

   template<typename T>
   __global__ void assigntable4d(const unsigned int i0, const unsigned int i1, const unsigned int d2, T* data, const unsigned int len)
   {
      newOne->assigntable4d(i0, i1, d2, data, len);
   }

   template<typename T>
   __global__ void assignx(const unsigned int i0, T* data, const unsigned int len)
   {
      newOne->assignx(i0, data, len);
   }

   template<typename T>
   __global__ void assigndomain(const unsigned int i0, T* data, const unsigned int len)
   {
      newOne->assigndomain(i0, data, len);
   }

   template<typename T>
   __global__ void assignrange(T* data, const unsigned int len)
   {
      newOne->assignrange(data, len);
   }

   __global__ void copydim(int d)
   {
      newOne->copydim(d);
   }

   __global__ void copytableFileName(const char* data)
   {
      newOne->copytableFileName(data);
   }

   __global__ void setInstance()
   {
      d_Factory->setInstance(newOne->gettableFileName(), newOne);
   }

   __global__ void resizetable2d(const unsigned int n)
   {
      newOne->resizetable2d(n);
   }

   __global__ void resizetable3d_0(const unsigned int d0)
   {
      newOne->resizetable3d_0(d0);
   }

   __global__ void resizetable3d_1(const unsigned int i0, const unsigned int d1)
   {
      newOne->resizetable3d_1(i0, d1);
   }

   __global__ void resizetable4d_0(const unsigned int d0)
   {
      newOne->resizetable4d_0(d0);
   }

   __global__ void resizetable4d_1(const unsigned int i0, const unsigned int d1)
   {
      newOne->resizetable4d_1(i0, d1);
   }

   __global__ void resizetable4d_2(const unsigned int i0, const unsigned int i1, const unsigned int d2)
   {
      newOne->resizetable4d_2(i0, i1, d2);
   }

   __global__ void resizex(const unsigned int n)
   {
      newOne->resizex(n);
   }

   __global__ void resizedomain(const unsigned int n)
   {
      newOne->resizedomain(n);
   }

   __global__ void allocOnCuda3di(float*** ptr, const unsigned int i, const unsigned int sz)
   {
      ptr[i] = new float*[sz];
   }

   __global__ void allocOnCuda3dij(float*** ptr, const unsigned int i, const unsigned int j, const unsigned int sz)
   {
      ptr[i][j] = new float[sz];
   }

   typedef float data_type;

   void transferDataToCuda(char const * tableFileName)
   {
      NUTableInterpolation const * ptr = getInstance(tableFileName);
      
      char* d_tableFileName = nullptr;
      checkCudaErrors(cudaMalloc((void **)&d_tableFileName, (ptr->gettableFileName().size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemset(d_tableFileName, NULL, (ptr->gettableFileName().size() + 1) * sizeof(char)));
      checkCudaErrors(cudaMemcpy(d_tableFileName, ptr->gettableFileName().c_str(), (ptr->gettableFileName().size() + 1) * sizeof(char), cudaMemcpyHostToDevice));
      initNewNUTableOnDevice << <1, 1 >> >(d_tableFileName);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      checkCudaErrors(cudaFree(d_tableFileName));

      data_type* d_table1d = nullptr;
      checkCudaErrors(cudaMalloc((void **)&d_table1d, ptr->gettable1d().size() * sizeof(data_type)));
      checkCudaErrors(cudaMemcpy(d_table1d, ptr->gettable1d().data(), ptr->gettable1d().size() * sizeof(data_type), cudaMemcpyHostToDevice));
      assigntable1d << <1, 1 >> >(d_table1d, ptr->gettable1d().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      data_type** d_table2d = new data_type*[ptr->gettable2d().size()];
      resizetable2d << <1, 1 >> >(ptr->gettable2d().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->gettable2d().size(); ++i) {
         checkCudaErrors(cudaMalloc((void **)&d_table2d[i], ptr->gettable2d()[i].size() * sizeof(data_type)));
         checkCudaErrors(cudaMemcpy(d_table2d[i], ptr->gettable2d()[i].data(), ptr->gettable2d()[i].size() * sizeof(data_type), cudaMemcpyHostToDevice));
         assigntable2d << <1, 1 >> >(i, d_table2d[i], ptr->gettable2d()[i].size());
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      delete[] d_table2d;

      data_type*** d_table3d = new data_type**[ptr->gettable3d().size()]; // 2d matrix of pointers on CPU each points to 1d vector of data on GPU
      resizetable3d_0 << <1, 1 >> >(ptr->gettable3d().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->gettable3d().size(); ++i) {
         d_table3d[i] = new data_type*[ptr->gettable3d()[i].size()];
         resizetable3d_1 << <1, 1 >> >(i, ptr->gettable3d()[i].size());
         for (int j = 0; j < ptr->gettable3d()[i].size(); ++j) {
            checkCudaErrors(cudaMalloc((void **)&d_table3d[i][j], ptr->gettable3d()[i][j].size() * sizeof(data_type)));
            checkCudaErrors(cudaMemcpy(d_table3d[i][j], ptr->gettable3d()[i][j].data(), ptr->gettable3d()[i][j].size() * sizeof(data_type), cudaMemcpyHostToDevice));
            assigntable3d << <1, 1 >> >(i, j, d_table3d[i][j], ptr->gettable3d()[i][j].size());
         }
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      for (int i = 0; i < ptr->gettable3d().size(); ++i) {
         delete[] d_table3d[i];
      }
      delete[] d_table3d;

      data_type**** d_table4d = new data_type***[ptr->gettable4d().size()];
      resizetable4d_0 << <1, 1 >> >(ptr->gettable4d().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->gettable4d().size(); ++i) {
         d_table4d[i] = new data_type**[ptr->gettable4d()[i].size()];
         resizetable4d_1 << <1, 1 >> >(i, ptr->gettable4d()[i].size());
         for (int j = 0; j < ptr->gettable4d()[i].size(); ++j) {
            d_table4d[i][j] = new data_type*[ptr->gettable4d()[i][j].size()];
            resizetable4d_2 << <1, 1 >> >(i, j, ptr->gettable4d()[i][j].size());
            for (int k = 0; k < ptr->gettable4d()[i][j].size(); ++k) {
               checkCudaErrors(cudaMalloc((void **)&d_table4d[i][j][k], ptr->gettable4d()[i][j][k].size() * sizeof(data_type)));
               checkCudaErrors(cudaMemcpy(d_table4d[i][j][k], ptr->gettable4d()[i][j][k].data(), ptr->gettable4d()[i][j][k].size() * sizeof(data_type), cudaMemcpyHostToDevice));
               assigntable4d << <1, 1 >> >(i, j, k, d_table4d[i][j][k], ptr->gettable4d()[i][j][k].size());
            }
         }
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->gettable4d().size(); ++i) {
         for (int j = 0; j < ptr->gettable4d()[i].size(); ++j) {
            delete[] d_table4d[i][j];
         }
         delete[] d_table4d[i];
      }
      delete[] d_table4d;

      data_type** d_x = new data_type*[ptr->getx().size()];
      resizex << <1, 1 >> >(ptr->getx().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->getx().size(); ++i) {
         checkCudaErrors(cudaMalloc((void **)&d_x[i], ptr->getx()[i].size() * sizeof(data_type)));
         checkCudaErrors(cudaMemcpy(d_x[i], ptr->getx()[i].data(), ptr->getx()[i].size() * sizeof(data_type), cudaMemcpyHostToDevice));
         assignx << <1, 1 >> >(i, d_x[i], ptr->getx()[i].size());
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      delete[] d_x;

      data_type** d_domain = new data_type*[ptr->getdomain().size()];
      resizedomain << <1, 1 >> >(ptr->getdomain().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      for (int i = 0; i < ptr->getdomain().size(); ++i) {
         checkCudaErrors(cudaMalloc((void **)&d_domain[i], ptr->getdomain()[i].size() * sizeof(data_type)));
         checkCudaErrors(cudaMemcpy(d_domain[i], ptr->getdomain()[i].data(), ptr->getdomain()[i].size() * sizeof(data_type), cudaMemcpyHostToDevice));
         assigndomain << <1, 1 >> >(i, d_domain[i], ptr->getdomain()[i].size());
      }
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
      delete[] d_domain;

      data_type* d_range = nullptr;
      checkCudaErrors(cudaMalloc((void **)&d_range, ptr->getrange().size() * sizeof(data_type)));
      checkCudaErrors(cudaMemcpy(d_range, ptr->getrange().data(), ptr->getrange().size() * sizeof(data_type), cudaMemcpyHostToDevice));
      assignrange << <1, 1 >> >(d_range, ptr->getrange().size());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      copydim << <1, 1 >> >(ptr->getdim());
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());

      setInstance << <1, 1 >> >();
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError());
   }
}
