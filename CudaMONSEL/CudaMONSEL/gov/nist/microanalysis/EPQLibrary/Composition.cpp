#include "Composition.h"

enum Composition::Representation
{
   UNDETERMINED, WEIGHT_PCT, STOICIOMETRY
};

const long Composition::serialVersionUID = 0x42;
const float Composition::OUT_OF_THIS_MANY_ATOMS = 1.0f;
