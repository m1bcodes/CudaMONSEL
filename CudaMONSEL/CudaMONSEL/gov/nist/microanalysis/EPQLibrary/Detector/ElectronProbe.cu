//import gov.nist.microanalysis.EPQLibrary.SpectrumProperties;
#include "gov\nist\microanalysis\EPQLibrary\ToSI.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace ElectronProbe
{
   //public ElectronProbe(SpectrumProperties props) {
   //   mProbeProperties = props.clone();
   //}

   //public ElectronProbe(String name) {
   //   super();
   //   mProbeProperties = new SpectrumProperties();
   //   mProbeProperties.setTextProperty(SpectrumProperties.Instrument, name);
   //}

   //public SpectrumProperties getProbeProperties() {
   //   return mProbeProperties;
   //}

   ///**
   //* Min beam energy in Joule
   //*
   //* @return double
   //*/
   //public double getMinBeamEnergy() {
   //   return mMinBeamEnergy;
   //}

   ///**
   //* Min beam energy in Joules
   //*/
   //public void setMinBeamEnergy(double minBeamEnergy) {
   //   mMinBeamEnergy = minBeamEnergy;
   //}

   ///**
   //* Max beam energy in Joule
   //*
   //* @return double
   //*/
   //public double getMaxBeamEnergy() {
   //   return mMaxBeamEnergy;
   //}

   ///**
   //* Max beam energy in Joules
   //*/
   //public void setMaxBeamEnergy(double maxBeamEnergy) {
   //   mMaxBeamEnergy = maxBeamEnergy;
   //}

   ///**
   //* Computes the position of the detector assuming that an axis from the
   //* optimal working distance (optWD) to the detector at an altitude angle
   //* (altitudeAngle) above the x-y plane for the specified distance will find
   //* the dispersive element of the detector. The detector is located at an
   //* orientation specified by azimuthAngle relative to the nominal x-axis.
   //*
   //* @param optWD
   //* @param altitudeAngle
   //* @param azimuthAngle
   //* @param distance
   //* @return double[3] with the detector position
   //*/
   void computePosition(const double optWD, const double altitudeAngle, const double azimuthAngle, const double distance, double res[])
   {
      res[0] = distance * ::cos(altitudeAngle) * ::cos(azimuthAngle);
      res[1] = distance * ::cos(altitudeAngle) * ::sin(azimuthAngle);
      res[2] = optWD - distance * ::sin(altitudeAngle);
   }

   ///**
   //* Computes the orientation of the detector assuming that the detector points
   //* directly towards the point (0,0,optWD) at an altitude and azimuth
   //* specified.
   //*
   //* @param optWD
   //* @param altitudeAngle
   //* @param azimuthAngle
   //* @param distance
   //* @return a double[3] normalized to 1
   //*/
   //public static double[] computeOrientation(double optWD, double altitudeAngle, double azimuthAngle, double distance) {
   //   return Math2.normalize(Math2.negative(computePosition(optWD, altitudeAngle, azimuthAngle, distance)));
   //}

   //@Override
   //   public String toString() {
   //   return mProbeProperties.getTextWithDefault(SpectrumProperties.Instrument, "Unknown");
   //}

   ///**
   //* Change the name of the ElectronProbe
   //*
   //* @param name
   //*/
   //public void setName(String name) {
   //   mProbeProperties.setTextProperty(SpectrumProperties.Instrument, name);
   //}

   //@Override
   //   public int compareTo(ElectronProbe o) {
   //   return toString().compareTo(o.toString());
   //}
}
