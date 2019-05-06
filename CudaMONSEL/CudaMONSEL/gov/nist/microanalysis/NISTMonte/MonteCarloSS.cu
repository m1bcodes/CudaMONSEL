#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\GaussianBeam.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"

#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace MonteCarloSS
{
   static const int ScatterEvent = 1;
   static const int NonScatterEvent = ScatterEvent + 1;
   static const int BackscatterEvent = ScatterEvent + 2;
   static const int ExitMaterialEvent = ScatterEvent + 3;
   static const int TrajectoryStartEvent = ScatterEvent + 4;
   static const int TrajectoryEndEvent = ScatterEvent + 5;
   static const int LastTrajectoryEvent = ScatterEvent + 6;
   static const int FirstTrajectoryEvent = ScatterEvent + 7;
   static const int StartSecondaryEvent = ScatterEvent + 8;
   static const int EndSecondaryEvent = ScatterEvent + 9;
   static const int PostScatterEvent = ScatterEvent + 10;

   static const int BeamEnergyChanged = 100;
   static const int XAxis = 0;
   static const int YAxis = 1;
   static const int ZAxis = 2;
   const float ChamberRadius = 0.1f;
   static const float SMALL_DISP = 1.0e-15f;

   MonteCarloSS::MonteCarloSS(ElectronGunT& gun, RegionT& chamber) : mGun(gun), mChamber(chamber)
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

   double distance(const double pos0[], const double pos1[])
   {
      return ::sqrt(Math2::sqr(pos1[0] - pos0[0]) + Math2::sqr(pos1[1] - pos0[1]) + Math2::sqr(pos1[2] - pos0[2]));
   }
}
