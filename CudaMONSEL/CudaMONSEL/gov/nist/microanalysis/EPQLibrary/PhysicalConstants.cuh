#ifndef _PHYSICAL_CONSTANTS_CUH_
#define _PHYSICAL_CONSTANTS_CUH_

#include <cuda_runtime.h>

namespace PhysicalConstants
{
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ extern const float AvagadroNumber; // dimensionless
   __constant__ extern const float SpeedOfLight; // m/s

   __constant__ extern const float ElectronCharge; // c
   __constant__ extern const float PlanckConstant; // j s
   __constant__ extern const float PlanckReduced; // j s
   __constant__ extern const float ElectronMass; // kg

   __constant__ extern const float ElectronRestMass; // j
   __constant__ extern const float ProtonMass; // kg
   __constant__ extern const float NeutronMass; // kg
   __constant__ extern const float UnifiedAtomicMass; // kg
   __constant__ extern const float PermittivityOfFreeSpace; // f/m
   __constant__ extern const float PermeabilityOfFreeSpace; // n/(a^2)
   __constant__ extern const float BoltzmannConstant; // j/k

   __constant__ extern const float GravitationConstant; // kg m^2
   __constant__ extern const float PlanckLength; // m
   __constant__ extern const float PlanckMass; // kg
   __constant__ extern const float PlanckTemperature; // dimensionless
   __constant__ extern const double PlanckTime; // s

   __constant__ extern const float RydbergEnergy; // joules

   __constant__ extern const float BohrRadius; // meter

   __constant__ extern const float FineStructure;

   __constant__ extern const float ClassicalElectronRadius;

   __constant__ extern const float IcePoint; // kelvin
   __constant__ extern const float StandardAtmosphere; // pascal

   __constant__ extern const float PI; // pi
#else
   extern const float AvagadroNumber; // dimensionless
   extern const float SpeedOfLight; // m/s

   extern const float ElectronCharge; // c
   extern const float PlanckConstant; // j s
   extern const float PlanckReduced; // j s
   extern const float ElectronMass; // kg

   extern const float ElectronRestMass; // j
   extern const float ProtonMass; // kg
   extern const float NeutronMass; // kg
   extern const float UnifiedAtomicMass; // kg
   extern const float PermittivityOfFreeSpace; // f/m
   extern const float PermeabilityOfFreeSpace; // n/(a^2)
   extern const float BoltzmannConstant; // j/k

   extern const float GravitationConstant; // kg m^2
   extern const float PlanckLength; // m
   extern const float PlanckMass; // kg
   extern const float PlanckTemperature; // dimensionless
   extern const double PlanckTime; // s

   extern const float RydbergEnergy; // joules

   extern const float BohrRadius; // meter

   extern const float FineStructure;

   extern const float ClassicalElectronRadius;

   extern const float IcePoint; // kelvin
   extern const float StandardAtmosphere; // pascal

   extern const float PI; // pi
#endif
};

#endif