#ifndef _MONTE_CARLO_SS_CU_
#define _MONTE_CARLO_SS_CU_

#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\NullMaterialScatterModel.cuh"
#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"
#include <vector>
#include <string>

////#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
////#include "gov\nist\microanalysis\NISTMonte\IMaterialScatterModel.cuh"
////#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
////#include "gov\nist\microanalysis\EPQLibrary\Element.cuh"
////#include "gov\nist\microanalysis\EPQLibrary\Material.cuh"
//
//#include <vector>
//#include <string>
//
//class MonteCarloSS
//{
//public:
//   typedef std::string MonteCarloSSNameT;
//
//   MonteCarloSS();
//
//   virtual int GetId() = 0;
//};

namespace MonteCarloSS
{
   extern const float ChamberRadius;

   class MonteCarloSS final
   {
   public:
      MonteCarloSS(ElectronGunT& gun, RegionT& mChamber);

   private:
      RegionT& mChamber;
      ElectronGunT& mGun;
      ElectronT * mElectron;
      std::vector<ElectronT*> mElectronStack;

      //virtual int GetId() = 0;
   };

   double distance(const double pos0[], const double pos1[]);

   static const NullMaterialScatterModelT NULL_MSM;
}

#endif
