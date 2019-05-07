// file: package gov\nist\microanalysis\EPQLibrary\Detector\ElectronProbe.cuh

#ifndef _ELECTRON_PROBE_CUH_
#define _ELECTRON_PROBE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace ElectronProbe
{
   // Properties that are shared by all detectors attached to the instrument
   //class ElectronProbe
   //{

   //private:
   //   SpectrumProperties mProbeProperties;
   //   double mMinBeamEnergy = ToSI.keV(5.0);
   //   double mMaxBeamEnergy = ToSI.keV(30.0);
   //};

   PositionVecT computePosition(double optWD, double altitudeAngle, double azimuthAngle, double distance);
}

#endif