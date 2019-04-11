#ifndef _PHYSICAL_CONSTANTS_CUH_
#define _PHYSICAL_CONSTANTS_CUH_

#include <cuda_runtime.h>

namespace PhysicalConstants
{
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
};

#endif