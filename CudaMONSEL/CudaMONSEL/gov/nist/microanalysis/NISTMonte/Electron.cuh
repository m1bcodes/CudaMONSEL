//#ifndef _ELECTRON_CUH_
//#define _ELECTRON_CUH_
//
//#include "gov\nist\microanalysis\NISTMonte\Types.cuh"
//
//namespace Electron
//{
//   class RegionBase;
//   class Element::Element;
//
//   class Electron
//   {
//   private:
//      // The x,y & z coordinates of the electron
//      PositionVecT mPosition; // final transient
//
//      // The location of the electron before the last call to updatePosition
//      PositionVecT mPrevPosition; // final transient
//
//      // The direction of the current trajectory segment
//      double mPhi, mTheta; // transient
//
//      // The kinetic energy of the electron
//      double mEnergy; // transient
//
//      // Kinetic energy of the electron upon conclusion of the previous step
//      double previousEnergy; // transient
//
//      int mStepCount; // transient
//
//      RegionBase& mCurrentRegion; // transient
//
//      RegionBase& mPrevRegion; // transient
//
//      Element::Element& mScatteringElement; // transient
//
//      bool mTrajectoryComplete; // transient
//
//      long ident; // A unique identifying number to assist tracking, final
//
//      long parentID = 0; // 0 if from e-gun. Otherwise ID of parent.
//   };
//
//   static long lastID = 0; // ID of last generated electron
//}
//
//#endif