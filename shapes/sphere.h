#ifndef SPHERE_H
#define SPHERE_H
#include "../vec3.h"
#include "../ray.h"
#include "../hit.h"
#include "../shape.h"
class Sphere : public Shape {
  public:
    Vec3 center;
    float radius;

    Sphere(const Vec3& _center, float _radius) : center(_center), radius(_radius) {};

    bool intersect(const Ray& ray, Hit& res) const {
      float a = ray.direction.length2();
      float b = 2*dot(ray.direction, ray.origin - center);
      float c = (ray.origin - center).length2() - radius*radius;
      float D = b*b - 4*a*c;
      if(D < 0) return false;

      float t0 = (-b - std::sqrt(D))/(2*a);
      float t1 = (-b + std::sqrt(D))/(2*a);

      float t = t0;
      if(t > ray.tmax) return false;
      if(t < ray.tmin) {
          t = t1;
          if(t < ray.tmin || t > ray.tmax) return false;
      }

      res.t = t;
      res.hitPos = ray(t);
      res.hitNormal = normalize(res.hitPos - center);
      res.hitSphere = this;

      float phi = std::atan2(res.hitNormal.z, res.hitNormal.x);
      if(phi < 0) phi += 2*M_PI;
      float theta = std::acos(res.hitNormal.y);
      res.u = phi/(2*M_PI);
      res.v = theta/M_PI;

      return true;
    };
};
#endif
