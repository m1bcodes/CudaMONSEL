#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\EPQLibrary\Detector\ElectronProbe.cuh"

#include "gov\nist\microanalysis\Utility\ActionListener.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\Random.cuh"

namespace MonteCarloSS
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const int ScatterEvent = 1;
   __constant__ const int NonScatterEvent = 2;
   __constant__ const int BackscatterEvent = 3;
   __constant__ const int ExitMaterialEvent = 4;
   __constant__ const int TrajectoryStartEvent = 5;
   __constant__ const int TrajectoryEndEvent = 6;
   __constant__ const int LastTrajectoryEvent = 7;
   __constant__ const int FirstTrajectoryEvent = 8;
   __constant__ const int StartSecondaryEvent = 9;
   __constant__ const int EndSecondaryEvent = 10;
   __constant__ const int PostScatterEvent = 11;

   __constant__ const int BeamEnergyChanged = 100;

   __constant__ static const int XAxis = 0;
   __constant__ static const int YAxis = 1;
   __constant__ static const int ZAxis = 2;
   __constant__ const float ChamberRadius = 0.1f;
   __constant__ const float SMALL_DISP = 1.0e-15f;
#else
   const int ScatterEvent = 1;
   const int NonScatterEvent = ScatterEvent + 1;
   const int BackscatterEvent = ScatterEvent + 2;
   const int ExitMaterialEvent = ScatterEvent + 3;
   const int TrajectoryStartEvent = ScatterEvent + 4;
   const int TrajectoryEndEvent = ScatterEvent + 5;
   const int LastTrajectoryEvent = ScatterEvent + 6;
   const int FirstTrajectoryEvent = ScatterEvent + 7;
   const int StartSecondaryEvent = ScatterEvent + 8;
   const int EndSecondaryEvent = ScatterEvent + 9;
   const int PostScatterEvent = ScatterEvent + 10;

   const int BeamEnergyChanged = 100;

   static const int XAxis = 0;
   static const int YAxis = 1;
   static const int ZAxis = 2;
   const float ChamberRadius = 0.1f;
   const float SMALL_DISP = 1.0e-15f;
#endif

   __host__ __device__ MonteCarloSS::MonteCarloSS(ElectronGunT const * gun, RegionT* chamber, ElectronT* electron) : mGun(gun), mChamber(chamber), mElectron(electron)
   {
      //TODO: shift the responsibility to the caller
      //const double center[] = {
      //   0.0,
      //   0.0,
      //   0.0
      //};
      //SphereT sphere(center, ChamberRadius);
      //mGun.setCenter(sphere.getInitialPoint().data());
      //mGun.setBeamEnergy(ToSI::keV(20.0));
      //mChamber = new RegionT(NULL, &NULL_MSM, &sphere);
   }

   //RegionT addSubRegion(const RegionT * parent, const MaterialT& mat, const ShapeT& shape)
   //{
   //   if (parent == NULL) printf("bad");
   //   return new RegionT(parent, new BasicMaterialModel(mat), shape);
   //}

   //RegionT* addSubRegion(RegionT& parent, IMaterialScatterModelT& msm, ShapeT& shape)
   //{
   //   return new RegionT(&parent, &msm, &shape); // TODO: deal with this, DO NOT USE IT
   //}

   //public Map<Material, Double> getMaterialMap(double[] startPt, double[] endPt) { // used
   //   final HashMap<Material, Double> traj = new HashMap<Material, Double>();
   //   double[] start = startPt;
   //   RegionBase region = mChamber.containingSubRegion(start);
   //   final double eps = 1.0e-7;
   //   while ((region != null) && (distance(start, endPt) > eps)) {
   //      final double[] end = endPt.clone();
   //      final RegionBase nextRegion = region.findEndOfStep(start, end);
   //      double dist = distance(start, end);
   //      if (dist > 0.0) {
   //         if (traj.containsKey(region.getMaterial()))
   //            dist += (traj.get(region.getMaterial())).doubleValue();
   //         traj.put(region.getMaterial(), new Double(dist));
   //      }
   //      start = Math2.plus(end, Math2.multiply(SMALL_DISP, Math2.normalize(Math2.minus(endPt, start))));
   //      region = nextRegion;
   //   }
   //   return traj;
   //}

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ static bool mDisableEvents = false;
#else
   static bool mDisableEvents = false;
#endif
   __host__ __device__ void MonteCarloSS::fireEvent(const int ae)
   {
      if (!(mEventListeners.empty() || mDisableEvents)) {
         for (auto sel : mEventListeners)
            sel->actionPerformed(ae);
      }
   }

   RegionT* MonteCarloSS::getChamber()
   {
      return mChamber;
   }

   __host__ __device__ void MonteCarloSS::initializeTrajectory()
   {
      mElectron = mGun->createElectron();
      auto reg = mChamber->containingSubRegion(mElectron->getPosition());
      mElectron->setCurrentRegion(reg);
      // Stop when you can't generate any more x-rays
      mElectron->setScatteringElement(nullptr);

   }

   __host__ __device__ void MonteCarloSS::takeStep()
   {
      const double *pos0 = mElectron->getPosition();

      const RegionBaseT* currentRegion = mElectron->getCurrentRegion();
      if ((currentRegion == nullptr) || !(currentRegion->getShape()->contains(pos0))) {
         currentRegion = mChamber->containingSubRegion(pos0);
         mElectron->setCurrentRegion(currentRegion);
         if (currentRegion == nullptr) {
            mElectron->setTrajectoryComplete(true);
            return;
         }
      }
      IMaterialScatterModelT* msm = currentRegion->getScatterModel();
      if (msm == nullptr) printf("MonteCarloSS::takeStep: msm is null\n");

      double pos1[3];
      mElectron->candidatePoint(msm->randomMeanPathLength(*mElectron), pos1);
      const RegionBaseT* nextRegion = currentRegion->findEndOfStep(pos0, pos1);
      mElectron->move(pos1, msm->calculateEnergyLoss(Math2::distance3d(pos0, pos1), *mElectron));
      const bool tc = (mElectron->getEnergy() < msm->getMinEforTracking()) || mElectron->isTrajectoryComplete();
      mElectron->setTrajectoryComplete(tc);
      if (!tc) {
         if (nextRegion == currentRegion) {
            if (mChamber == nullptr) printf("MonteCarloSS::takeStep: mChamber == nullptr");
            if (mElectron == nullptr) printf("MonteCarloSS::takeStep: mElectron == nullptr");
            if (currentRegion == nullptr) printf("MonteCarloSS::takeStep: currentRegion == nullptr");
            fireEvent(ScatterEvent);
            ElectronT* secondary = msm->scatter(*mElectron);
            fireEvent(PostScatterEvent);
            mElectron->setTrajectoryComplete((mElectron->getEnergy() < msm->getMinEforTracking()) || mElectron->isTrajectoryComplete());
            if (secondary != nullptr) trackSecondaryElectron(secondary);

            if (mElectron->getCurrentRegion() != currentRegion) printf("MonteCarloSS::takeStep: mElectron->getCurrentRegion() != currentRegion\n");
         }
         else if (nextRegion != nullptr) {
            fireEvent(NonScatterEvent);
            ElectronT* secondary = msm->barrierScatter(mElectron, nextRegion);
            double candpt[3];
            mElectron->candidatePoint(SMALL_DISP, candpt);
            mElectron->setPosition(candpt);
            if (!(mElectron->getCurrentRegion()->getShape()->contains(mElectron->getPosition())))
               mElectron->setCurrentRegion(mChamber->containingSubRegion(mElectron->getPosition()));
            if (mElectron->getCurrentRegion() != currentRegion)
               fireEvent(ExitMaterialEvent);
            if (secondary != nullptr) {
               double secandpt[3];
               secondary->candidatePoint(SMALL_DISP, secandpt);
               secondary->setPosition(secandpt);
               trackSecondaryElectron(secondary);
            }
         }
         else {
            fireEvent(BackscatterEvent);
            mElectron->setCurrentRegion(nullptr);
            mElectron->setTrajectoryComplete(true);
         }
      }
   }

   __host__ __device__ void MonteCarloSS::trackSecondaryElectron(ElectronT* newElectron)
   {
      double mMinEnergy = newElectron->getCurrentRegion()->getScatterModel()->getMinEforTracking();
      if (newElectron->getEnergy() > mMinEnergy) {
         fireEvent(StartSecondaryEvent);
         mElectronStack.push(mElectron);
         mElectron = newElectron;
         fireEvent(StartSecondaryEvent);
      }
   }

   //int MonteCarloSS::getElectronGeneration() const
   //{
   //   return mElectronStack.size();
   //}

   __host__ __device__ void MonteCarloSS::addActionListener(ActionListenerT& sel)
   {
      mEventListeners.push_back(&sel);
   }

   __host__ __device__ void MonteCarloSS::removeActionListener(ActionListenerT& sel)
   {
      auto itr = amp::find(mEventListeners.begin(), mEventListeners.end(), &sel);
      if (itr != mEventListeners.end()) {
         mEventListeners.erase(itr);
      }
   }

   __host__ __device__ bool MonteCarloSS::allElectronsComplete()
   {
      bool tc = mElectron->isTrajectoryComplete();
      while (tc && !mElectronStack.empty()) {
         fireEvent(EndSecondaryEvent);
         delete mElectron;
         mElectron = mElectronStack.top();
         mElectronStack.pop();

         tc = mElectron->isTrajectoryComplete();
      }
      if (tc && mElectronStack.empty() && mElectron) delete mElectron; // deleting the last electron
      return tc;
   }

   __host__ __device__ void MonteCarloSS::runTrajectory()
   {
      initializeTrajectory();
      fireEvent(TrajectoryStartEvent);
      //unsigned int i = 0;
      while (!allElectronsComplete()) {
         takeStep();
         //++i;
      }
      fireEvent(TrajectoryEndEvent);
      //printf("%d steps\n", i);
   }

   __host__ __device__ void MonteCarloSS::runMultipleTrajectories(int n)
   {
      fireEvent(FirstTrajectoryEvent);
      for (int i = 0; i < n; ++i) {
         //printf("itr #%d:\n", i);
         runTrajectory();
      }
      fireEvent(LastTrajectoryEvent);
   }

   __host__ __device__ double MonteCarloSS::getBeamEnergy() const
   {
      return mGun->getBeamEnergy();
   }

   void MonteCarloSS::setElectronGun(ElectronGunT& gun)
   {
      gun.setBeamEnergy(mGun->getBeamEnergy());
      mGun = &gun;
   }

   //void MonteCarloSS::setBeamEnergy(double beamEnergy)
   //{
   //   //mGun->setBeamEnergy(beamEnergy);
   //   //fireEvent(BeamEnergyChanged);
   //}

   void MonteCarloSS::computeDetectorPosition(double elevation, double theta, double res[])
   {
      double frac = 0.999;
      double r = frac * ChamberRadius;
      //if (mChamber->mShape instanceof Sphere)
      //   r = frac * ((Sphere)mChamber.mShape).getRadius();
      ElectronProbe::computePosition(0.0, elevation, theta, r, res);
   }

   __host__ __device__ const ElectronT& MonteCarloSS::getElectron() const
   {
      return *mElectron;
   }

   //public Set<AtomicShell> getAtomicShellSet()
   //   throws EPQException{ // not used
   //   final Set<AtomicShell> res = new TreeSet<AtomicShell>();
   //   final Set<Element> elements = mChamber.getElements(true);
   //   for (final Object element : elements) {
   //      final Element el = (Element)element;
   //      for (int sh = AtomicShell.K; sh < AtomicShell.NI; ++sh) {
   //         final AtomicShell shell = new AtomicShell(el, sh);
   //         final double ee = shell.getEdgeEnergy();
   //         if ((ee > 0.0) && (ee < mGun.getBeamEnergy()) && (XRayTransition.getStrongestLine(shell) != null))
   //            res.add(shell);
   //      }
   //   }
   //   return res;
   //}

   Element::UnorderedSetT MonteCarloSS::getElementSet() const
   {
      return mChamber->getElements(true);
   }

   //public void estimateTrajectoryVolume(double[] c0, double[] c1) {
   //   c0[0] = (c0[1] = (c0[2] = Double.MAX_VALUE));
   //   c1[0] = (c1[1] = (c1[2] = -Double.MAX_VALUE));
   //   final int nTraj = 100;
   //   mDisableEvents = true;
   //   for (int i = 0; i < nTraj; ++i) {
   //      initializeTrajectory();
   //      while (!mElectron.isTrajectoryComplete()) {
   //         takeStep();
   //         final double[] endPt = mElectron.getPosition();
   //         final RegionBase endRegion = mChamber.containingSubRegion(endPt);
   //         if ((endRegion != null) && (endRegion != mChamber))
   //            for (int j = 0; j < 3; ++j) {
   //               if (endPt[j] < c0[j])
   //                  c0[j] = endPt[j];
   //               if (endPt[j] > c1[j])
   //                  c1[j] = endPt[j];
   //            }
   //      }
   //   }
   //   mDisableEvents = false;
   //}

   //void updateMaterial(const MaterialT& oldMat, const MaterialT& newMat)
   //{
   //   mChamber.updateMaterial(oldMat, new BasicMaterialModel(newMat));
   //}

   void MonteCarloSS::rotate(double pivot[], double phi, double theta, double psi)
   {
      for (const RegionBaseT* r : mChamber->getSubRegions()) {
         //if (r instanceof TransformableRegion) ((TransformableRegionT*)r)->rotate(pivot, phi, theta, psi);
      }
   }

   void MonteCarloSS::translate(double distance[])
   {
      for (const RegionBaseT* r : mChamber->getSubRegions()) {
         //if (r instanceof TransformableRegion) ((TransformableRegionT*)r)->translate(distance);
      }
   }

   const RegionBaseT* MonteCarloSS::findRegionContaining(double point[]) const
   {
      return mChamber->containingSubRegion(point);
   }
}
