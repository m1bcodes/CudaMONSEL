#ifndef _BACKSCATTER_STATS_CUH_
#define _BACKSCATTER_STATS_CUH_

#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\ActionListener.cuh"

namespace BackscatterStats
{
   class Datum
   {
   public:
      Datum(long eID, long tStep, double e0, const double pos[], double theta, double phi);
      Datum(const Datum&);

      long getElectronID() const;
      long getTrajStep() const;
      double getEnergyeV() const;
      VectorXd getPosition() const;
      double getTheta() const;
      double getPhi() const;
      StringT toString() const;

      static StringT getHeader();

   private:
      const long electronID;
      const long trajStep;
      const double mkEeV;
      const VectorXd mPosition;
      const double mTheta;
      const double mPhi;
   };

   class BackscatterStats : public ActionListenerT
   {
   public:
      BackscatterStats(const MonteCarloSST& mcss);
      BackscatterStats(const MonteCarloSST& mcss, int nEnergyBins);
      ~BackscatterStats();

      void actionPerformed(const int ae) override;

      const HistogramT& backscatterEnergyHistogram() const;
      const HistogramT& forwardscatterEnergyHistogram() const;
      const HistogramT& elevationHistogram() const;
      const HistogramT& azimuthalHistogram() const;

      double backscatterFraction() const;
      int getEnergyBinCount() const;
      double BackscatterStats::forwardscatterFraction() const;
      double BackscatterStats::scatterFraction() const;
      bool BackscatterStats::getLogDetected() const;
      void BackscatterStats::setLogDetected(bool logDetected);
      std::vector<Datum> getLog() const;

   private:
      void initialize();

      int mEnergyBinCount;
      const MonteCarloSST& mMonte;
      double mBeamEnergy; // in eV
      HistogramT* mFwdEnergyBins;
      HistogramT* mBackEnergyBins;
      HistogramT* mAzimuthalBins;
      HistogramT* mElevationBins;
      int mEventCount;

      bool mLogDetected;
      std::vector<Datum> mLog;
   };
}

#endif