#include "rasterizer.h"
#include <random> 
using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    int sqrt_rate = (int)std::sqrt(sample_rate);  
    for (int i = 0; i < sqrt_rate; ++i) {
      for (int j = 0; j < sqrt_rate; ++j) {
        size_t sample_index = ((size_t)y * width + x) * sample_rate + (j * sqrt_rate + i);
        sample_buffer[sample_index] = c;
      }
    }
}
  
  
  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  void RasterizerImp::rasterize_triangle(float x0, float y0,
                                        float x1, float y1,
                                        float x2, float y2,
                                        Color color) {
    // Vertex data arrays
    float px[3] = { x0, x1, x2 };
    float py[3] = { y0, y1, y2 };

    // Compute signed area of triangle
    float total = (px[1] - px[0]) * (py[2] - py[0]) -
                  (px[2] - px[0]) * (py[1] - py[0]);
    if (fabs(total) < 1e-8f) return;  // Degenerate triangle

    // Make orientation Clockwise
    if (total < 0) {
      std::swap(px[1], px[2]);
      std::swap(py[1], py[2]);
      total = -total;
    }

    // Bounding box in screen pixels
    int x_lo = std::max(0, (int)floor(std::min({ px[0], px[1], px[2] })));
    int x_hi = std::min((int)width, (int)ceil(std::max({ px[0], px[1], px[2] })));
    int y_lo = std::max(0, (int)floor(std::min({ py[0], py[1], py[2] })));
    int y_hi = std::min((int)height, (int)ceil(std::max({ py[0], py[1], py[2] })));

    // Supersampling configuration
    int ss = (int)sqrt(sample_rate);
    float delta = 1.0f / ss;

    for (int y = y_lo; y < y_hi; ++y) {
      for (int x = x_lo; x < x_hi; ++x) {
        for (int s = 0; s < sample_rate; ++s) {
          int sub_x = s % ss;
          int sub_y = s / ss;

          float sx = x + (sub_x + 0.5f) * delta;
          float sy = y + (sub_y + 0.5f) * delta;

          // Compute edge-function area signs
          float a0 = (px[1] - sx) * (py[2] - sy) - (px[2] - sx) * (py[1] - sy);
          float a1 = (px[2] - sx) * (py[0] - sy) - (px[0] - sx) * (py[2] - sy);
          float a2 = (px[0] - sx) * (py[1] - sy) - (px[1] - sx) * (py[0] - sy);
          if (a0 >= 0 && a1 >= 0 && a2 >= 0) {
            size_t idx = ((size_t)y * width + x) * sample_rate + s;
            sample_buffer[idx] = color;
          }
        }
      }
    }
  }


  // Computes barycentric coordinates (alpha, beta, gamma) for point (px, py)
  // with respect to triangle defined by vertices (x0, y0), (x1, y1), (x2, y2)
  void compute_barycentric_coordinates(
    float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    float px, float py,
    float &alpha, float &beta, float &gamma) {

    // Compute vectors
    float v0x = x1 - x0;
    float v0y = y1 - y0;
    float v1x = x2 - x0;
    float v1y = y2 - y0;
    float v2x = px - x0;
    float v2y = py - y0;

    // Compute dot products
    float d00 = v0x * v0x + v0y * v0y;
    float d01 = v0x * v1x + v0y * v1y;
    float d11 = v1x * v1x + v1y * v1y;
    float d20 = v2x * v0x + v2y * v0y;
    float d21 = v2x * v1x + v2y * v1y;

    // Compute denominator
    float denom = d00 * d11 - d01 * d01;

    // Compute barycentric coordinates
    beta = (d11 * d20 - d01 * d21) / denom;
    gamma = (d00 * d21 - d01 * d20) / denom;
    alpha = 1.0f - beta - gamma;
  }

  void RasterizerImp::rasterize_interpolated_color_triangle(
    float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2) {

    Vector3D z(0, 0, 1);
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);

    if (cross(((p1 + p2) / 2) - p0, p1 - p0).z < 0) {
      swap(p1, p2);
      swap(c1, c2);
    }

    Vector3D lin0 = p0 - p1;
    Vector3D lin1 = p1 - p2;
    Vector3D lin2 = p2 - p0;
    Vector3D n0 = cross(z, lin0);
    Vector3D n1 = cross(z, lin1);
    Vector3D n2 = cross(z, lin2);

    // Bounding box limits for triangle
  float x_list[] = { static_cast<float>(p0.x), static_cast<float>(p1.x), static_cast<float>(p2.x) };
  float y_list[] = { static_cast<float>(p0.y), static_cast<float>(p1.y), static_cast<float>(p2.y) };

    // x pixel bounds
    int minx = std::max((int)floor(*std::min_element(x_list, x_list + 3)), 0);
    int maxx = std::min((int)ceil(*std::max_element(x_list, x_list + 3)), (int)width);
    // y pixel bounds
    int miny = std::max((int)floor(*std::min_element(y_list, y_list + 3)), 0);
    int maxy = std::min((int)ceil(*std::max_element(y_list, y_list + 3)), (int)height);

    int sr = sqrt(sample_rate);   


    Matrix3x3 M(p0.x, p1.x, p2.x, p0.y, p1.y, p2.y, 1, 1, 1);
    M = M.inv();

    for (int y = miny; y < maxy; y++) {
      for (int x = minx; x < maxx; x++) {
        int s = 0;
        for (int j = 0; j < sr; j++) {
          float py = (float)y + ((float)j + 0.5f) / (float)sr;
          for (int i = 0; i < sr; i++) {
            float px = (float)x + ((float)i + 0.5f) / (float)sr;
            Vector3D p(px, py, 1);

            if ((dot(p - p1, n0) >= 0) && (dot(p - p2, n1) >= 0) && (dot(p - p0, n2) >= 0)) {
              Vector3D weights = M * p;
              weights.z = 1 - weights.x - weights.y;
              
              sample_buffer[sample_rate * (y * width + x) + s] = weights.x * c0 + weights.y * c1 + weights.z * c2;
            }
            s++;
          }
        }
      }
    }
  }

  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
                                                float x1, float y1, float u1, float v1,
                                                float x2, float y2, float u2, float v2,
                                                Texture& tex) {
    Vector3D A(x0, y0, 0), B(x1, y1, 0), C(x2, y2, 0), Z(0, 0, 1);
    if (cross((B + C) * 0.5 - A, B - A).z < 0) std::swap(B, C);

    Vector3D N0 = cross(Z, A - B);
    Vector3D N1 = cross(Z, B - C);
    Vector3D N2 = cross(Z, C - A);

    float x_vals[3] = { x0, x1, x2 }, y_vals[3] = { y0, y1, y2 };
    int x_start = floor(std::min({ x0, x1, x2 }));
    int x_end   = ceil(std::max({ x0, x1, x2 }));
    int y_start = floor(std::min({ y0, y1, y2 }));
    int y_end   = ceil(std::max({ y0, y1, y2 }));

    int res = (int)std::sqrt(sample_rate);
    float step = 1.0f / res;

    Matrix3x3 Bmat(x0, x1, x2, y0, y1, y2, 1, 1, 1);
    Bmat = Bmat.inv();
    Vector3D U(u0, u1, u2), V(v0, v1, v2);

    SampleParams samp;
    samp.psm = psm;
    samp.lsm = lsm;

    for (int y = y_start; y < y_end; ++y) {
        for (int x = x_start; x < x_end; ++x) {
            int idx = 0;
            for (int j = 0; j < res; ++j) {
                float sy = y + (j + 0.5f) * step;
                for (int i = 0; i < res; ++i) {
                    float sx = x + (i + 0.5f) * step;
                    Vector3D P(sx, sy, 1);
                    if (dot(P - B, N0) >= 0 && dot(P - C, N1) >= 0 && dot(P - A, N2) >= 0) {
                        Vector3D W = Bmat * P;
                        Vector3D W_dx = Bmat * (P + Vector3D(1, 0, 0));
                        Vector3D W_dy = Bmat * (P + Vector3D(0, 1, 0));

                        W.z = 1.0f - W.x - W.y;
                        W_dx.z = 1.0f - W_dx.x - W_dx.y;
                        W_dy.z = 1.0f - W_dy.x - W_dy.y;

                        samp.p_uv = Vector2D(dot(W, U), dot(W, V));
                        samp.p_dx_uv = Vector2D(dot(W_dx, U), dot(W_dx, V)) - samp.p_uv;
                        samp.p_dy_uv = Vector2D(dot(W_dy, U), dot(W_dy, V)) - samp.p_uv;

                        Color col = tex.sample(samp);
                        fill_supersample(x, y, idx, col);
                    }
                    ++idx;
                }
            }
        }
    }
  }



  void RasterizerImp::fill_supersample(size_t px,        // pixel-x
    size_t py,        // pixel-y
    size_t sub,       // sub-sample index
    Color  col)       // colour to deposit
  {
  /* fast early–out if any coordinate is outside the buffer */
  if (px >= width || py >= height || sub >= sample_rate)
  return;

  /* flatten (px,py,sub) -> linear index and store the colour */
  size_t offset = ((py * width + px) * sample_rate) + sub;
  sample_buffer[offset] = col;
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    this->sample_rate = rate;
    sample_buffer.resize(width * height * sample_rate, Color::White);
}
  

  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height) {
    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
}

  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        Color color(0, 0, 0);
        for (int i = 0; i < sample_rate; ++i) {
          color += sample_buffer[(y * width + x) * sample_rate + i];
        }
        color.r /= sample_rate;
        color.g /= sample_rate;
        color.b /= sample_rate;
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&color.r)[k] * 255;
        }
      }
    }
  } 
  
  Rasterizer::~Rasterizer() { }


}// CGL
