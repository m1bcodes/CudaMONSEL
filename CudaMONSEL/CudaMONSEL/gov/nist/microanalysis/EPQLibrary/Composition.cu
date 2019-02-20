#include "Composition.cuh"

namespace Composition
{
   __device__ const long serialVersionUID = 0x42;
   __device__ const double OUT_OF_THIS_MANY_ATOMS = 1.0;

   enum Representation
   {
      UNDETERMINED, WEIGHT_PCT, STOICIOMETRY
   };
}
