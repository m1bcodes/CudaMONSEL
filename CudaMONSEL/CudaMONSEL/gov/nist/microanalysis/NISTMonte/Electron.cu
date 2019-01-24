#ifndef ELECTRON_CU
#define ELECTRON_CU

#include "cuda_runtime.h"
//#include "MonteCarloSS.cu"

class Electron
{
private:
   float mPosition[3], mPrevPosition[3];
   float mPhi, mTheta;
   float mEnergy, mPrevEnergy;
   int mStepCount;

   //MonteCarloSS::RegionBase mCurrentRegion, mPrevRegion;
   //Element mScatteringElement;

   int id;

public:
   __host__ __device__ Electron(int id) { this->id = id; }
   __host__ __device__ ~Electron() { };

   __host__ __device__ int GetId() { return id; }
};

#endif
