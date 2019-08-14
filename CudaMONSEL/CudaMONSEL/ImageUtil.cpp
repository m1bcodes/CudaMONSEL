#include "ImageUtil.h"

#include <windows.h>
#include <fstream>

#include <gdiplus.h>

namespace ImageUtil
{
   bool saveImage(const char* szPathName, const char* lpBits, int w, int h)
   {
      HDC hDC = GetDC(NULL);
      HDC memHDC = CreateCompatibleDC(hDC);
      BITMAPINFO* pbmi = (BITMAPINFO*) new BYTE[sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD) * 256];
      memset(pbmi, 0, sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD) * 256);
      BITMAPINFO& bmi = *pbmi;
      bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
      bmi.bmiHeader.biWidth = w;
      bmi.bmiHeader.biHeight = -h; // top-down
      bmi.bmiHeader.biPlanes = 1;
      bmi.bmiHeader.biBitCount = 8;
      bmi.bmiHeader.biCompression = BI_RGB;
      bmi.bmiHeader.biSizeImage = (((w * bmi.bmiHeader.biBitCount + 31) & ~31) >> 3) * h;
      bmi.bmiHeader.biClrUsed = 256; // size of palette, 256 for grayscale

      for (int i = 0; i < 256; i++) {
         bmi.bmiColors[i].rgbRed = i;
         bmi.bmiColors[i].rgbGreen = i;
         bmi.bmiColors[i].rgbBlue = i;
         bmi.bmiColors[i].rgbReserved = 0;
      }
      HBITMAP bitmap = CreateDIBSection(hDC, &bmi, DIB_RGB_COLORS, (void**)&lpBits[0], NULL, NULL);
      ReleaseDC(NULL, hDC);
      DeleteDC(hDC);

      SelectObject(memHDC, bitmap);

      BITMAPFILEHEADER bf;
      memset(&bf, 0, sizeof(BITMAPFILEHEADER));
      bf.bfType = MAKEWORD('B', 'M');
      bf.bfOffBits = sizeof(BITMAPFILEHEADER) + bmi.bmiHeader.biSize + sizeof(RGBQUAD) * 256;
      //bf.bfOffBits = sizeof(BITMAPFILEHEADER) + bmi.bmiHeader.biSize;
      bf.bfSize = bf.bfOffBits + bmi.bmiHeader.biSizeImage;

      std::vector<char> bitmapDataV;
      bitmapDataV.insert(bitmapDataV.end(), (char*)&bf, ((char*)&bf) + sizeof(BITMAPFILEHEADER));
      bitmapDataV.insert(bitmapDataV.end(), (char*)&bmi.bmiHeader, ((char*)&bmi.bmiHeader) + sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD) * 256);
      bitmapDataV.insert(bitmapDataV.end(), &lpBits[0], &lpBits[0] + bmi.bmiHeader.biSizeImage);

      DeleteObject(SelectObject(memHDC, bitmap));
      DeleteObject(bitmap);
      ReleaseDC(NULL, memHDC);
      DeleteDC(memHDC);

      std::ofstream pFile(szPathName, std::ios_base::binary);
      if (!pFile.is_open()) {
         return false;
      }

      pFile.write(&bitmapDataV[0], bitmapDataV.size());
      pFile.close();
      delete[] (BYTE*)pbmi;

      return true;
   }

   void saveResults(const std::string path, const unsigned int* data, int w, int h)
   {
      std::vector<char> bitmapData(w * h, 0);
      for (int k = 0; k < h; ++k) {
         for (int l = 0; l < w; ++l) {
            bitmapData[k * w + l] = data[k * w + l];
         }
      }

      saveImage(path.c_str(), bitmapData.data(), w, h);
   }

   void drawGradient()
   {
      int img_x = 256;
      int img_y = 256;
      std::vector<char> bitmapData(img_x * img_y, 0);
      for (int k = 0; k < img_y; ++k) {
         for (int l = 0; l < img_x; ++l) {
            bitmapData[k * img_x + l] = k * img_x + l;
            std::cout << (int)bitmapData[k * img_x + l] << " ";
         }
         std::cout << std::endl;
      }

      saveImage("f.bmp", bitmapData.data(), img_x, img_y);
   }
}