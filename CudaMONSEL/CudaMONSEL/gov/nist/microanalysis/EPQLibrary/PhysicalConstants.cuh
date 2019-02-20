#ifndef _PHYSICAL_CONSTANTS_CUH_
#define _PHYSICAL_CONSTANTS_CUH_

#include <cuda_runtime.h>

namespace PhysicalConstants
{
   extern __device__ const double AvagadroNumber; // dimensionless
   extern __device__ const double SpeedOfLight; // m/s

   extern __device__ const double ElectronCharge; // c
   extern __device__ const double PlanckConstant; // j s
   extern __device__ const double PlanckReduced; // j s
   extern __device__ const double ElectronMass; // kg

   extern __device__ const double ElectronRestMass; // j
   extern __device__ const double ProtonMass; // kg
   extern __device__ const double NeutronMass; // kg
   extern __device__ const double UnifiedAtomicMass; // kg
   extern __device__ const double PermittivityOfFreeSpace; // f/m
   extern __device__ const double PermeabilityOfFreeSpace; // n/(a^2)
   extern __device__ const double BoltzmannConstant; // j/k

   extern __device__ const double GravitationConstant; // kg m^2
   extern __device__ const double PlanckLength; // m
   extern __device__ const double PlanckMass; // kg
   extern __device__ const double PlanckTemperature; // dimensionless
   extern __device__ const double PlanckTime; // s

   extern __device__ const double RydbergEnergy; // joules

   extern __device__ const double BohrRadius; // meter

   extern __device__ const double FineStructure;

   extern __device__ const double ClassicalElectronRadius;

   extern __device__ const double IcePoint; // kelvin
   extern __device__ const double StandardAtmosphere; // pascal
};

#endif