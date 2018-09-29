#ifndef SPHERE_H
#define SPHERE_H
#include "../vec3.h"
#include "../ray.h"
#include "../hit.h"
#include "../shape.h"
class Sphere : public Shape {
  public:
    Vec3 center;
    double radius;

    Sphere(const std::shared_ptr<Material>& _material, const std::string& _type, const Vec3& _center, double _radius) : Shape(_material, _type), center(_center), radius(_radius) {};

    bool intersect(const Ray& ray, Hit& res) const {
      double a = ray.direction.length2();
      double b = 2*dot(ray.direction, ray.origin - center);
      double c = (ray.origin - center).length2() - radius*radius;
      double D = b*b - 4*a*c;
      if(D < 0) return false;

      double t0 = (-b - std::sqrt(D))/(2*a);
      double t1 = (-b + std::sqrt(D))/(2*a);

      double t = t0;
      if(t > ray.tmax) return false;
      if(t < ray.tmin) {
        t = t1;
        if(t < ray.tmin || t > ray.tmax) return false;
      }

      res.t = t;
      res.hitPos = ray(t);
      res.hitNormal = normalize(res.hitPos - center);
      res.hitShape = this;

      double phi = std::atan2(res.hitNormal.z, res.hitNormal.x);
      if(phi < 0) phi += 2*M_PI;
      double theta = std::acos(res.hitNormal.y);
      res.u = phi/(2*M_PI);
      res.v = theta/M_PI;

      return true;
    };
};
#endif
