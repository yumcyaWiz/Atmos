#ifndef RAY_H
#define RAY_H
#include "vec3.h"
class Ray {
  public:
    Vec3 origin;
    Vec3 direction;
    mutable double tmin;
    mutable double tmax;

    Ray() : tmin(0.001f), tmax(1e9) {};
    Ray(const Vec3& origin, const Vec3& direction) : tmin(0.001f), tmax(1e9), origin(origin), direction(direction) {};

    Vec3 operator()(double t) const {
      return origin + t*direction;
    };
};
#endif
