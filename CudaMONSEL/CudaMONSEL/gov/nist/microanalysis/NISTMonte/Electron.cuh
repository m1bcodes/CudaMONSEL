#ifndef _ELECTRON_CUH_
#define _ELECTRON_CUH_

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
   Electron(int id) { this->id = id; }
   //~Electron() { };

   int GetId() { return id; }
};

#endif
