#include "gov\nist\microanalysis\NISTMonte\BackscatterStats.cuh"

#include "gov\nist\microanalysis\EPQLibrary\FromSI.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\Utility\Histogram.cuh"

namespace BackscatterStats
{
   long Datum::getElectronID() const
   {
      return electronID;
   }

   long Datum::getTrajStep() const
   {
      return trajStep;
   }

   double Datum::getEnergyeV() const
   {
      return mkEeV;
   }

   const double* Datum::getPosition() const
   {
      return mPosition;
   }

   double Datum::getTheta() const
   {
      return mTheta;
   }

   double Datum::getPhi() const
   {
      return mPhi;
   }

   __host__ __device__ Datum::Datum(long eID, long tStep, double e0, const double pos[], double theta, double phi) :
      electronID(eID),
      trajStep(tStep),
      mkEeV(e0),
      mTheta(theta),
      mPhi(phi)
   {
      memcpy(mPosition, pos, sizeof(mPosition[0]) * 3);
   }

   __host__ __device__ Datum::Datum() :
      electronID(0),
      trajStep(0),
      mkEeV(0),
      mTheta(0),
      mPhi(0)
   {
      mPosition[0] = 0;
      mPosition[1] = 0;
      mPosition[2] = 0;
   }

   __host__ __device__ Datum::Datum(const Datum& other) :
      electronID(other.electronID),
      trajStep(other.trajStep),
      mkEeV(other.mkEeV),
      mTheta(other.mTheta),
      mPhi(other.mPhi)
   {
      memcpy(mPosition, other.mPosition, sizeof(mPosition[0]) * 3);
   }

   __host__ __device__ Datum& Datum::operator=(const Datum& other)
   {
      electronID = (other.electronID);
      trajStep = (other.trajStep);
      mkEeV = (other.mkEeV);
      mTheta = (other.mTheta);
      mPhi = (other.mPhi);
      memcpy(mPosition, other.mPosition, sizeof(mPosition[0]) * 3);

      return *this;
   }

   __host__ __device__ bool Datum::operator==(const Datum& other)
   {
      if (this == &other) return true;
      return
         electronID == (other.electronID) &&
         trajStep == (other.trajStep) &&
         mkEeV == (other.mkEeV) &&
         mTheta == (other.mTheta) &&
         mPhi == (other.mPhi) &&
         mPosition[0] == mPosition[0] &&
         mPosition[1] == mPosition[1] &&
         mPosition[2] == mPosition[2];
   }

   StringT Datum::toString() const
   {
      StringT sb = 
         "BackscatterStats: " +
         amp::to_string(electronID) + "\t" + 
         amp::to_string(trajStep) + "\t" +
         amp::to_string(mkEeV) + "\t" +
         amp::to_string(mPosition[0]) + "\t" +
         amp::to_string(mPosition[1]) + "\t" +
         amp::to_string(mPosition[2]) + "\t" +
         amp::to_string(mTheta) + "\t" +
         amp::to_string(mPhi);
      //sb.append(electronID);
      //sb.append("\t");
      //sb.append(trajStep);
      //sb.append("\t");
      //sb.append(mkEeV);
      //sb.append("\t");
      //sb.append(mPosition[0]);
      //sb.append("\t");
      //sb.append(mPosition[1]);
      //sb.append("\t");
      //sb.append(mPosition[2]);
      //sb.append("\t");
      //sb.append(mTheta);
      //sb.append("\t");
      //sb.append(mPhi);
      return sb;
   }

   StringT Datum::getHeader()
   {
      return "Electron ID\tTraj Step\tkinetic E (eV)\tx\ty\tz\ttheta\tphi";
   }

   BackscatterStats::BackscatterStats(const MonteCarloSST& mcss) :
      mEnergyBinCount(400),
      mMonte(mcss),
      mBeamEnergy(FromSI::eV(mMonte.getBeamEnergy())),
      mEventCount(0),
      mFwdEnergyBins(new HistogramT(0.0, mBeamEnergy, mEnergyBinCount)),
      mBackEnergyBins(new HistogramT(0.0, mBeamEnergy, mEnergyBinCount)),
      mAzimuthalBins(new HistogramT(0.0, 2.0 * Math2::PI, 360)),
      mElevationBins(new HistogramT(0.0, Math2::PI, 180)),
      mLogDetected(false)
   {
   }

   __host__ __device__ BackscatterStats::BackscatterStats(const MonteCarloSST& mcss, int nEnergyBins) :
      mEnergyBinCount(nEnergyBins),
      mMonte(mcss),
      mBeamEnergy(FromSI::eV(mMonte.getBeamEnergy())),
      mFwdEnergyBins(new HistogramT(0.0, mBeamEnergy, mEnergyBinCount)),
      mBackEnergyBins(new HistogramT(0.0, mBeamEnergy, mEnergyBinCount)),
      mAzimuthalBins(new HistogramT(0.0, 2.0 * Math2::PI, 360)),
      mElevationBins(new HistogramT(0.0, Math2::PI, 180)),
      mEventCount(0),
      mLogDetected(false)
   {
   }

   __host__ __device__ BackscatterStats::~BackscatterStats()
   {
      delete mFwdEnergyBins;
      delete mBackEnergyBins;
      delete mAzimuthalBins;
      delete mElevationBins;
   }

   __host__ __device__ void BackscatterStats::initialize()
   {
      mBeamEnergy = FromSI::eV(mMonte.getBeamEnergy());
      mFwdEnergyBins->clear();
      mBackEnergyBins->clear();
      mAzimuthalBins->clear();
      mElevationBins->clear();
      //mFwdEnergyBins = new HistogramT(0.0, mBeamEnergy, mEnergyBinCount);
      //mBackEnergyBins = new HistogramT(0.0, mBeamEnergy, mEnergyBinCount);
      mLog.clear();
   }

   __host__ __device__ void BackscatterStats::actionPerformed(const int ae)
   {
      //assert(ae.getSource() instanceof MonteCarloSS);
      //assert(ae.getSource() == mMonte);
      if (ae == MonteCarloSS::FirstTrajectoryEvent) {
         mEventCount = 0;
      }
      else if (ae == MonteCarloSS::BackscatterEvent) {
         const ElectronT& el = mMonte.getElectron();
         //auto pos = el.getPosition();
         const double* pos = el.getPosition();
         const double elevation = (Math2::PI / 2) - ::atan2(pos[2], ::sqrt((pos[0] * pos[0]) + (pos[1] * pos[1])));
         if (!(elevation >= 0.0)) printf("BackscatterStats::actionPerformed: !(elevation >= 0.0) %.10e\n", elevation);
         if (!(elevation <= Math2::PI)) printf("BackscatterStats::actionPerformed: !(elevation <= Math2::PI) %.10e\n", elevation);
         double azimuth = ::atan2(pos[1], pos[0]);
         if (azimuth < 0.0)
            azimuth = (2.0 * Math2::PI) + azimuth;
         if (!(azimuth >= 0.0)) printf("BackscatterStats::actionPerformed: !(azimuth >= 0.0) %.10e\n", azimuth);
         if (!(azimuth <= 2.0 * Math2::PI)) printf("BackscatterStats::actionPerformed: !(azimuth >= 0.0) %.10e\n", azimuth);
         //synchronized(this) {
         mElevationBins->add(elevation);
         mAzimuthalBins->add(azimuth);
         double kEeV = FromSI::eV(el.getEnergy());
         if (kEeV > FromSI::eV(el.getInitialEnergy())) {
            printf("Histogram::add: energy higher than energy of eletron spawned (%.5e/%.5e)\n", kEeV, FromSI::eV(el.getInitialEnergy()));
            kEeV = FromSI::eV(el.getInitialEnergy());
         }
         if (elevation < (Math2::PI / 2.0))
            mFwdEnergyBins->add(kEeV);
         else
            mBackEnergyBins->add(kEeV);
         if (mLogDetected) {
            mLog.push_back(Datum(el.getIdent(), el.getStepCount(), kEeV, pos, el.getTheta(), el.getPhi()));
         }
         //}
      }
      else if (ae == MonteCarloSS::TrajectoryEndEvent) {
         //synchronized(this) {
         mEventCount++;
         //}
      }
      else if (ae == MonteCarloSS::BeamEnergyChanged) {
         //synchronized(this) {
         initialize();
         //}
      }
   }

   __host__ __device__ const HistogramT& BackscatterStats::backscatterEnergyHistogram() const
   {
      return *mBackEnergyBins;
   }

   __host__ __device__ const HistogramT& BackscatterStats::forwardscatterEnergyHistogram() const
   {
      return *mFwdEnergyBins;
   }

   const HistogramT& BackscatterStats::elevationHistogram() const
   {
      return *mElevationBins;
   }

   const HistogramT& BackscatterStats::azimuthalHistogram() const
   {
      return *mAzimuthalBins;
   }

   //void dump(final OutputStream os)
   //{
   //   final NumberFormat nf = NumberFormat.getInstance();
   //   nf.setMaximumFractionDigits(3);
   //   final PrintWriter pw = new PrintWriter(os);
   //   { // Header
   //      pw.println("Beam energy\t" + nf.format(mBeamEnergy / 1000.0) + " keV");
   //      pw.println("Backscatter\t" + Integer.toString(mBackEnergyBins.totalCounts()));
   //      pw.println("Forwardscatter\t" + Integer.toString(mFwdEnergyBins.totalCounts()));
   //   }
   //   pw.println("Forward and back scattered electron energy histogram");
   //   pw.println("Bin\tBack\tForward");
   //   assert mBackEnergyBins.binCount() == mFwdEnergyBins.binCount();
   //   final Iterator<Integer> bs = mBackEnergyBins.getResultMap("{0,number,#.##}").values().iterator();
   //   for (final Map.Entry<Histogram.BinName, Integer> me : mFwdEnergyBins.getResultMap("{0,number,#.##}").entrySet()) {
   //      pw.print(me.getKey().toString());
   //      pw.print("\t");
   //      pw.print(bs.next().toString());
   //      pw.print("\t");
   //      pw.println(me.getValue().toString());
   //   }

   //   pw.println("Azimuthal angle histogram");
   //   pw.println("Bin\tAngle");
   //   for (final Map.Entry<Histogram.BinName, Integer> me : mAzimuthalBins.getResultMap("{0,number,#.##}").entrySet()) {
   //      pw.print(me.getKey().toString());
   //      pw.print("\t");
   //      pw.println(me.getValue().toString());
   //   }

   //   pw.println("Elevation angle histogram");
   //   pw.println("Bin\tAngle");
   //   for (final Map.Entry<Histogram.BinName, Integer> me : mElevationBins.getResultMap("{0,number,#.##}").entrySet()) {
   //      pw.print(me.getKey().toString());
   //      pw.print("\t");
   //      pw.println(me.getValue().toString());
   //   }

   //   /* If logging is turned on, output data for each detected electron */
   //   if (mLogDetected) {
   //      pw.println("Detected electron log (electron ID, energy, position, and direction of motion at detection)");
   //      pw.println("Number of logged electrons: " + Integer.toString(mLog.size()));
   //      pw.println(Datum.getHeader());
   //      for (final Datum logEntry : mLog)
   //         pw.println(logEntry.toString());
   //   }
   //   pw.close();
   //}

   __host__ __device__ double BackscatterStats::backscatterFraction() const
   {
      return (double)mBackEnergyBins->totalCounts() / (double)mEventCount;
   }

   int BackscatterStats::getEnergyBinCount() const
   {
      return mEnergyBinCount;
   }

   double BackscatterStats::forwardscatterFraction() const
   {
      return (double)mFwdEnergyBins->totalCounts() / (double)mEventCount;
   }

   double BackscatterStats::scatterFraction() const
   {
      return backscatterFraction() + forwardscatterFraction();
   }

   bool BackscatterStats::getLogDetected() const
   {
      return mLogDetected;
   }

   void BackscatterStats::setLogDetected(bool logDetected)
   {
      mLogDetected = logDetected;
   }

   const amp::vector<Datum>& BackscatterStats::getLog() const
   {
      return mLog;
   }
}
