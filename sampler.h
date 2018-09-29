#ifndef SAMPLER_H
#define SAMPLER_H
#include "vec2.h"
class Sampler {
  public:
    virtual double getNext() = 0;
    virtual Vec2 getNext2D() = 0;
};


inline Vec3 sampleSphere(const Vec2& u) {
  const double z = 1 - 2*u.x;
  const double r = std::sqrt(std::max(0.0, 1.0 - z*z));
  const double phi = 2*M_PI*u.y;
  return Vec3(r*std::cos(phi), z, r*std::sin(phi));
}


#endif
