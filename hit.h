#ifndef HIT_H
#define HIT_H
#include "vec3.h"
#include "ray.h"

class Sphere;

class Hit {
  public:
    float t;
    Vec3 hitPos;
    Vec3 hitNormal;
    const Sphere* hitSphere;
    float u;
    float v;

    Hit() {};
};
#endif
