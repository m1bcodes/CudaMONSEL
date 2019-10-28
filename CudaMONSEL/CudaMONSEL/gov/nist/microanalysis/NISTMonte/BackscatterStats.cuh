#ifndef _BACKSCATTER_STATS_CUH_
#define _BACKSCATTER_STATS_CUH_

#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\ActionListener.cuh"

namespace BackscatterStats
{
   class Datum
   {
   public:
      __host__ __device__ Datum(long eID, long tStep, double e0, const double pos[], double theta, double phi);
      __host__ __device__ Datum(); // TODO: remove this (ie update vector class so that it does not rely on objects having a default constructor)
      __host__ __device__ Datum(const Datum&);

      __host__ __device__ Datum& operator=(const Datum&); // TODO: remove this (ie update vector class)
      __host__ __device__ bool operator==(const Datum&);

      long getElectronID() const;
      long getTrajStep() const;
      double getEnergyeV() const;
      const double* getPosition() const;
      double getTheta() const;
      double getPhi() const;
      StringT toString() const;

      StringT getHeader();

   private:
      long electronID;
      long trajStep;
      double mkEeV;
      double mPosition[3];
      double mTheta;
      double mPhi;
   };

   class BackscatterStats : public ActionListenerT
   {
   public:
      BackscatterStats(const MonteCarloSST& mcss);
      __host__ __device__ BackscatterStats(const MonteCarloSST& mcss, int nEnergyBins);
      __host__ __device__ ~BackscatterStats();

      __host__ __device__ void actionPerformed(const int ae) override;

      __host__ __device__ const HistogramT& backscatterEnergyHistogram() const;
      __host__ __device__ const HistogramT& forwardscatterEnergyHistogram() const;
      const HistogramT& elevationHistogram() const;
      const HistogramT& azimuthalHistogram() const;

      __host__ __device__ double backscatterFraction() const;
      int getEnergyBinCount() const;
      double forwardscatterFraction() const;
      double scatterFraction() const;
      bool getLogDetected() const;
      void setLogDetected(bool logDetected);
      const amp::vector<Datum>& getLog() const;

   private:
      __host__ __device__ void initialize();

      int mEnergyBinCount;
      const MonteCarloSST& mMonte;
      double mBeamEnergy; // in eV
      HistogramT* mFwdEnergyBins;
      HistogramT* mBackEnergyBins;
      HistogramT* mAzimuthalBins;
      HistogramT* mElevationBins;
      int mEventCount;

      bool mLogDetected;
      amp::vector<Datum> mLog;
   };
}

#endif