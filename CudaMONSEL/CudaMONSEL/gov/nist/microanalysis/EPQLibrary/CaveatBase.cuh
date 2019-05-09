// file: gov\nist\microanalysis\EPQLibrary\CaveatBase.cuh
#ifndef _CAVEAT_BASE_CUH_
#define _CAVEAT_BASE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace CaveatBase
{
   extern const char None[];
   extern const char Broken[];
   extern const char NotImplemented[];

   StringT append(StringT base, StringT str);
}

#endif