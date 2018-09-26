#ifndef SAMPLER_H
#define SAMPLER_H
#include "vec2.h"
class Sampler {
  public:
    virtual float getNext() = 0;
    virtual Vec2 getNext2D() = 0;
};


inline Vec3 sampleSphere(const Vec2& u) {
  const float z = 1 - 2*u.x;
  const float r = std::sqrt(std::max(0.0f, 1.0f - z*z));
  const float phi = 2*M_PI*u.y;
  return Vec3(r*std::cos(phi), z, r*std::sin(phi));
}


#endif
