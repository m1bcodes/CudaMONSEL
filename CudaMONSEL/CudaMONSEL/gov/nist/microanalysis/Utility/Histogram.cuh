// file: gov\nist\microanalysis\Utility\Histogram.cuh

#ifndef _HISTOGRAM_CUH_
#define _HISTOGRAM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Histogram
{
   //class BinName
   //{
   //public:
   //   __host__ __device__ BinName(int, const char[], const HistogramT&);

   //   StringT toString() const;

   //   __host__ __device__ bool operator<(const BinName& o) const;

   //private:
   //   const int mBin;
   //   const StringT mFormat;
   //   const HistogramT& mEnclosingClass;
   //};

   //typedef amp::unordered_map<BinName*, int> BinMap;

   class Histogram
   {
   public:
      __host__ __device__ Histogram(double min, double max, int nBins);
      __host__ __device__ Histogram(double min, double max, double ratio);
      __host__ __device__ Histogram(double binMins[], int len, double max);
      __host__ __device__ Histogram(const Histogram& hist);

      __host__ __device__ void addBin(double binMin);
      __host__ __device__ int bin(double val) const;
      double minValue(int bin) const;
      double maxValue(int bin) const;
      __host__ __device__ void add(double val);
      void add(double vals[], int len);
      __host__ __device__ int binCount() const;
      StringT binName(int bin, StringT format) const;
      __host__ __device__ int counts(int bin) const;
      __host__ __device__ void clear();
      int overrange() const;
      int underrange() const;
      __host__ __device__ int totalCounts() const;
      bool isBinMin(double binMin) const;
      void removeBin(int binNum);

      //BinMap getResultMap(const char format[]);

   private:
      VectorXd mBinMin;
      VectorXi mCounts;
   };
}

#endif