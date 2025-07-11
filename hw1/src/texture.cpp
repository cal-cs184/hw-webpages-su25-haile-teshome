#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    if (mipmap.empty()) return Color(1, 0, 1);
    int max_level = (int)mipmap.size() - 1;
    float level = 0;
    if (sp.lsm == L_ZERO) {
        level = 0;
    } else {
        level = get_level(sp);
        if (sp.lsm == L_NEAREST) level = round(level);
    }
    if (level < 0) level = 0;
    if (level > max_level) level = max_level;
    if (sp.lsm == L_ZERO || sp.lsm == L_NEAREST) {
        int ilevel = (int)level;
        if (sp.psm == P_NEAREST) {
            return sample_nearest(sp.p_uv, ilevel);
        } else {
            return sample_bilinear(sp.p_uv, ilevel);
        }
    } else if (sp.lsm == L_LINEAR) {
        int l0 = (int)floor(level);
        int l1 = std::min(l0 + 1, max_level);
        float w = level - l0;
        Color c0, c1;
        if (sp.psm == P_NEAREST) {
            c0 = sample_nearest(sp.p_uv, l0);
            c1 = sample_nearest(sp.p_uv, l1);
        } else {
            c0 = sample_bilinear(sp.p_uv, l0);
            c1 = sample_bilinear(sp.p_uv, l1);
        }
        return (1 - w) * c0 + w * c1;
    }
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    Vector2D du_dx = (sp.p_dx_uv - sp.p_uv) * (float)mipmap[0].width;
    Vector2D du_dy = (sp.p_dy_uv - sp.p_uv) * (float)mipmap[0].height;
    float Lx = sqrt(du_dx.x * du_dx.x + du_dx.y * du_dx.y);
    float Ly = sqrt(du_dy.x * du_dy.x + du_dy.y * du_dy.y);
    float L = std::max(Lx, Ly);
    float level = log2f(L);
    if (level < 0) level = 0;
    if (level > (int)mipmap.size() - 1) level = (int)mipmap.size() - 1;
    return level;
  }

  Color MipLevel::get_texel(int tx, int ty)const {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    auto& mip = mipmap[level];
    int tx = static_cast<int>(round(uv.x * (mip.width - 1)));
    int ty = static_cast<int>(round(uv.y * (mip.height - 1)));
      if (tx < 0 || tx >= mip.width || ty < 0 || ty >= mip.height) {
      return Color(1, 0, 1);  
    }
  
    return mip.get_texel(tx, ty);
  }
  
  Color Texture::sample_bilinear(Vector2D uv, int level) {
    auto& mip = mipmap[level];
    float tx = uv.x * (mip.width - 1);
    float ty = uv.y * (mip.height - 1);
  
    int x0 = static_cast<int>(floor(tx));
    int x1 = x0 + 1;
    int y0 = static_cast<int>(floor(ty));
    int y1 = y0 + 1;
  
    float sx = tx - x0;
    float sy = ty - y0;
  
    if (x0 < 0 || y0 < 0 || x1 >= mip.width || y1 >= mip.height)
      return Color(1, 0, 1);
  
    Color c00 = mip.get_texel(x0, y0);
    Color c10 = mip.get_texel(x1, y0);
    Color c01 = mip.get_texel(x0, y1);
    Color c11 = mip.get_texel(x1, y1);
  
    // Horizontal interpolation
    Color h0 = c00 * (1.0f - sx) + c10 * sx;
    Color h1 = c01 * (1.0f - sx) + c11 * sx;
  
    // Vertical interpolation
    Color result = h0 * (1.0f - sy) + h1 * sy;
    return result;
  }
  
  Color lerp(const Color& a, const Color& b, float t) {
    return a * (1.0f - t) + b * t;
}



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
