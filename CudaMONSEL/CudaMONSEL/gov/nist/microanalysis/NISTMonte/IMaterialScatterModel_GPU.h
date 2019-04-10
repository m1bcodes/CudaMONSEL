//#ifndef I_MATERIAL_SCATTER_MODEL_H
//#define I_MATERIAL_SCATTER_MODEL_H
//
////#include "../EPQLibrary/Material.cu"
////#include "Electron.cu"
//
//class IMaterialScatterModel
//{
//public:
//   virtual Material getMaterial() = 0;
//   virtual double getMinEforTracking() = 0;
//   virtual void setMinEforTracking(double minEforTracking) = 0;
//   virtual double randomMeanPathLength(Electron pe) = 0;
//   virtual Electron scatter(Electron pe) = 0;
//   virtual Electron barrierScatter(Electron pe, RegionBase nextRegion) = 0;
//   virtual double calculateEnergyLoss(double len, Electron pe) = 0;
//}
//
//#endif