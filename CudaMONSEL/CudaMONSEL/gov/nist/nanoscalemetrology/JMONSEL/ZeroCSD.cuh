// file: gov\nist\nanoscalemetrology\JMONSEL\ZeroCSD.cuh

#include "gov\nist\nanoscalemetrology\JMONSEL\SlowingDownAlg.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"

namespace ZeroCSD
{
   class ZeroCSD : public SlowingDownAlgT
   {
   public:
      ZeroCSD();

      double compute(double d, const ElectronT* pe) const override;
      void setMaterial(const SEmaterialT* mat) override;

      StringT toString() const override;
   };
}