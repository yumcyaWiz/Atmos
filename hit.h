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
    int iteration;

    Hit() : t(1e9), iteration(0) {};
};
#endif
