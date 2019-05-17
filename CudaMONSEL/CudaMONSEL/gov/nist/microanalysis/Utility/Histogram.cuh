// file: gov\nist\microanalysis\Utility\Histogram.cuh

#ifndef _HISTOGRAM_CUH_
#define _HISTOGRAM_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace Histogram
{
   class BinName
   {
   public:
      BinName(int, const char[], const HistogramT&);

      StringT toString() const;

      bool operator<(const BinName& o) const;

   private:
      const int mBin;
      const StringT mFormat;
      const HistogramT& mEnclosingClass;
   };

   typedef std::map<BinName*, int> BinMap;

   class Histogram
   {
   public:
      Histogram(double min, double max, int nBins);
      Histogram(double min, double max, double ratio);
      Histogram(double binMins[], int len, double max);
      Histogram(const Histogram& hist);

      void addBin(double binMin);
      int bin(double val) const;
      double minValue(int bin) const;
      double maxValue(int bin) const;
      void add(double val);
      void add(double vals[], int len);
      int binCount() const;
      StringT binName(int bin, StringT format) const;
      int counts(int bin) const;
      void clear();
      int overrange() const;
      int underrange() const;
      int totalCounts() const;
      bool isBinMin(double binMin) const;
      void removeBin(int binNum);

      //BinMap getResultMap(const char format[]);

   private:
      VectorXd mBinMin;
      VectorXi mCounts;
   };
}

#endif