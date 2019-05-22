#ifndef _NSHAPE_CUH_
#define _NSHAPE_CUH_

#include "gov\nist\microanalysis\NISTMonte\Declarations.cuh"

namespace NShapes
{
   extern NormalShapeT* createLine(
      double topz, // z of the top face
      double width, // line width
      double length, // length of line
      double thetal, // angle of right sidewall
      double thetar, // angle of left sidewall
      double radl, // radius of top right corner
      double radr // radius of top left corner
      );
}

#endif