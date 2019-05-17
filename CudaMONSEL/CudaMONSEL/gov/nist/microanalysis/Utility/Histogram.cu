#include "gov\nist\microanalysis\Utility\Histogram.cuh"

#include "Amphibian\Algorithm.cuh"

namespace Histogram
{
   BinName::BinName(int bin, const char format[], const Histogram& ec) : mBin(bin), mFormat(format), mEnclosingClass(ec)
   {
      if (!(bin >= -1)) printf("BinName::BinName: !(bin >= -1)\n");
      if (!(bin < mEnclosingClass.binCount() + 1)) printf("BinName::BinName: !(bin < binCount() + 1)\n");
   }

   StringT BinName::toString() const
   {
      StringT first = mBin < 0 ? "-inf" : mFormat + std::to_string(mEnclosingClass.minValue(mBin));
      StringT second = mBin >= mEnclosingClass.binCount() ? "inf" : mFormat + std::to_string(mEnclosingClass.maxValue(mBin));
      return "[" + first + "-" + second  + ")";
   }

   bool BinName::operator<(const BinName& rhs) const
   {
      return mBin < rhs.mBin;
   }
   
   Histogram::Histogram(double min, double max, int nBins) : mCounts(nBins + 2), mBinMin(nBins + 1)
   {
      if (max < min) {
         const double tmp = min;
         min = max;
         max = tmp;
      }
      if (min == max)
         printf("Histogram: min can not equal max\n");
      if (!(min < max)) printf("Histogram::Histogram: !(min < max)\n");
      if (!(nBins > 0)) printf("Histogram::Histogram: !(nBins > 0)\n");
      
      const double delta = (max - min) / nBins;
      mBinMin[0] = min;
      for (int i = 1; i < mBinMin.size(); ++i)
         mBinMin[i] = min + i * delta;
      Algorithm::quicksort(mBinMin.data(), 0, mBinMin.size() - 1);
   }

   Histogram::Histogram(double min, double max, double ratio) : mCounts((int)(::log(max / min) / ::log(ratio)) + 2), mBinMin((int)(::log(max / min) / ::log(ratio)) + 1)
   {
      if (ratio <= 1.0) printf("Histogram: ration must be greater than 1.0\n");
      if (min >= max) printf("Histogram: min must be less than max\n");
      if (min <= 0.0) printf("Histogram: min must be larger than 0.0\n");
      for (int i = 0; i < mBinMin.size(); ++i, min *= ratio) mBinMin[i] = min;
   }

   Histogram::Histogram(double binMins[], int len, double max) : mBinMin(len + 1), mCounts(len + 2)
   {
      for (int i = 0; i < len; ++i) mBinMin[i] = binMins[i];
      mBinMin[len] = max;
      Algorithm::quicksort(mBinMin.data(), 0, mBinMin.size() - 1);
      if (mBinMin[len - 1] != max) printf("Histogram: Max is not larger than all binMins.");
   }

   Histogram::Histogram(const Histogram& hist) : mBinMin(hist.mBinMin), mCounts(hist.mCounts)
   {
   }

   void Histogram::addBin(double binMin)
   {
      VectorXd newBinMin(mBinMin.begin(), mBinMin.end());
      newBinMin.push_back(binMin);
      Algorithm::quicksort(newBinMin.data(), 0, newBinMin.size() - 1);
      mCounts.resize(newBinMin.size() + 2);
      mBinMin = newBinMin;
   }

   int Histogram::bin(double val) const
   {
      int i = Algorithm::binarySearch(mBinMin.data(), 0, mBinMin.size() - 1, val);
      i = (i >= 0 ? i : -i - 2);
      if (!(i >= -1)) printf("Histogram::bin: index is %s for %s\n", std::to_string(i).c_str(), std::to_string(val).c_str());
      if (!(i < mCounts.size())) printf("Histogram::bin: !(i < mCounts.size())\n");
      return i;
   }

   double Histogram::minValue(int bin) const
   {
      return bin > -1 ? mBinMin[bin] : -INFINITY;
   }

   double Histogram::maxValue(int bin) const
   {
      return bin + 1 < mBinMin.size() ? mBinMin[bin + 1] : INFINITY;
   }

   void Histogram::add(double val)
   {
      ++mCounts[bin(val) + 1];
   }

   void Histogram::add(double vals[], int len)
   {
      for (int i = 0; i < len; ++i)
         add(vals[i]);
   }

   int Histogram::binCount() const
   {
      return mBinMin.size() - 1;
   }

   StringT Histogram::binName(int bin, StringT format) const
   {
      return (BinName(bin, format.c_str(), *this)).toString();
   }

   int Histogram::counts(int bin) const
   {
      return mCounts[bin + 1];
   }

   void Histogram::clear()
   {
      mCounts.assign(mCounts.size(), 0);
   }

   int Histogram::overrange() const
   {
      return mCounts[mCounts.size() - 1];
   }

   int Histogram::underrange() const
   {
      return mCounts[0];
   }

   int Histogram::totalCounts() const
   {
      int res = 0;
      for (auto c : mCounts)
         res += c;
      return res;
   }

   bool Histogram::isBinMin(double binMin) const
   {
      const int i = Algorithm::binarySearch(mBinMin.data(), 0, mBinMin.size() - 1, binMin);
      return i >= 0;
   }

   void Histogram::removeBin(int binNum)
   {
      if ((binNum < 0) || (binNum > mCounts.size() - 2)) printf("Histogram::removeBin: %d is not a bin in the histogram\n", binNum);
      VectorXd newBinMin(mBinMin.size() - 1);
      for (int index = 0; index < mBinMin.size(); index++)
         if (index < binNum)
            newBinMin[index] = mBinMin[index];
         else if (index > binNum)
            newBinMin[index - 1] = mBinMin[index];
      mBinMin = newBinMin;
      Algorithm::quicksort(mBinMin.data(), 0, mBinMin.size() - 1);

      VectorXi newCounts(mCounts.size() - 1);
      mCounts[binNum] += mCounts[binNum + 1];
      for (int index = 0; index < mCounts.size(); index++)
         if (index < binNum + 1)
            newCounts[index] = mCounts[index];
         else if (index > binNum + 1)
            newCounts[index - 1] = mCounts[index];
      mCounts = newCounts;
   }

   //BinMap Histogram::getResultMap(const char format[])
   //{
   //   BinMap res;
   //   for (int i = -1; i < binCount() + 1; ++i)
   //      res[new BinName(i, format, *this)] = counts(i); // TODO: fix
   //   return res;
   //}
}
