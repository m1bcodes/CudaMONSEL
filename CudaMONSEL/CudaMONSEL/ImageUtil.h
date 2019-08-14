#ifndef _IMAGE_UTIL_H_
#define _IMAGE_UTIL_H_

#include <iostream>
#include <vector>

namespace ImageUtil
{
   bool saveImage(const char* szPathName, const char* lpBits, int w, int h);
   void saveResults(const std::string path, const unsigned int* data, int w, int h);
   void drawGradient();
}

#endif