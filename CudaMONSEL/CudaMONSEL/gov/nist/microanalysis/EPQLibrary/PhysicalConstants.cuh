#ifndef _PHYSICAL_CONSTANTS_CUH_
#define _PHYSICAL_CONSTANTS_CUH_

#include <cuda_runtime.h>

namespace PhysicalConstants
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ extern const double AvagadroNumber; // dimensionless
   __constant__ extern const double SpeedOfLight; // m/s

   __constant__ extern const double ElectronCharge; // c
   __constant__ extern const double PlanckConstant; // j s
   __constant__ extern const double PlanckReduced; // j s
   __constant__ extern const double ElectronMass; // kg

   __constant__ extern const double ElectronRestMass; // j
   __constant__ extern const double ProtonMass; // kg
   __constant__ extern const double NeutronMass; // kg
   __constant__ extern const double UnifiedAtomicMass; // kg
   __constant__ extern const double PermittivityOfFreeSpace; // f/m
   __constant__ extern const double PermeabilityOfFreeSpace; // n/(a^2)
   __constant__ extern const double BoltzmannConstant; // j/k

   __constant__ extern const double GravitationConstant; // kg m^2
   __constant__ extern const double PlanckLength; // m
   __constant__ extern const double PlanckMass; // kg
   __constant__ extern const double PlanckTemperature; // dimensionless
   __constant__ extern const double PlanckTime; // s

   __constant__ extern const double RydbergEnergy; // joules

   __constant__ extern const double BohrRadius; // meter

   __constant__ extern const double FineStructure;

   __constant__ extern const double ClassicalElectronRadius;

   __constant__ extern const double IcePoint; // kelvin
   __constant__ extern const double StandardAtmosphere; // pascal

   __constant__ extern const double PI; // pi
#else
   extern const double AvagadroNumber; // dimensionless
   extern const double SpeedOfLight; // m/s

   extern const double ElectronCharge; // c
   extern const double PlanckConstant; // j s
   extern const double PlanckReduced; // j s
   extern const double ElectronMass; // kg

   extern const double ElectronRestMass; // j
   extern const double ProtonMass; // kg
   extern const double NeutronMass; // kg
   extern const double UnifiedAtomicMass; // kg
   extern const double PermittivityOfFreeSpace; // f/m
   extern const double PermeabilityOfFreeSpace; // n/(a^2)
   extern const double BoltzmannConstant; // j/k

   extern const double GravitationConstant; // kg m^2
   extern const double PlanckLength; // m
   extern const double PlanckMass; // kg
   extern const double PlanckTemperature; // dimensionless
   extern const double PlanckTime; // s

   extern const double RydbergEnergy; // joules

   extern const double BohrRadius; // meter

   extern const double FineStructure;

   extern const double ClassicalElectronRadius;

   extern const double IcePoint; // kelvin
   extern const double StandardAtmosphere; // pascal

   extern const double PI; // pi
#endif
};

#endif