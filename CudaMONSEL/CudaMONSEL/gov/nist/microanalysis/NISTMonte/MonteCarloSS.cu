#ifndef MONTE_CARLO_SS_CU
#define MONTE_CARLO_SS_CU

#include "cuda_runtime.h"

#include "Electron.cu"

class MonteCarloSS
{
public:
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
   static const float ChamberRadius;
   static const float SMALL_DISP;

public:
   class ElectronGun;
   class Shape;
   class RegionBase;
   class TransformableRegion;

   __host__ __device__ MonteCarloSS() { }
   __host__ __device__ ~MonteCarloSS() { };

   __host__ __device__ virtual int GetId() = 0;
};

#endif