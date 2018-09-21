#ifndef TERRAIN_H
#define TERRAIN_H
#include <cmath>
#include <memory>
#include "../vec2.h"
#include "../shape.h"
#include "sphere.h"


class Terrain : public Shape {
  private:
    std::shared_ptr<Sphere> sphere;

    float glslRandom(const Vec2& u) const {
      float i;
      return std::modf(std::sin(u.x * 12.9898 + u.y * 78.233) * 43758.5453123, &i);
    };
    float heightMap(const Vec2& u) const {
      return (std::sin(100*u.x)*std::sin(100*u.y) + 1);
    };
    Vec3 heightMapNormal(const Vec3& p) const {
      Vec3 n = normalize(p - sphere->center);
      Vec2 u = approxUV(p);

      float phi = 2*M_PI*u.x;
      float theta = M_PI*u.y;
      float px = sphere->radius * std::cos(phi) * std::sin(theta);
      float py = sphere->radius * std::cos(theta);
      float pz = sphere->radius * std::sin(phi) * std::sin(theta);
      Vec3 dpdu = Vec3(-2*M_PI*pz, 0, 2*M_PI*px);
      Vec3 dpdv = M_PI * Vec3(py*std::cos(phi), -sphere->radius*std::sin(theta), py*std::sin(phi));
      return normalize(cross(dpdu, dpdv));

      float du = (heightMap(u + Vec2(0.00001, 0)) - heightMap(u - Vec2(0.00001, 0)));
      float dv = (heightMap(u + Vec2(0, 0.00001)) - heightMap(u - Vec2(0, 0.00001)));
      Vec3 grad = normalize(du * dpdu + dv * dpdv);
      return grad;
    };

    Vec2 approxUV(const Vec3& p) const {
      Vec3 dir = normalize(p - sphere->center);
      float phi = std::atan2(dir.z, dir.x);
      if(phi < 0) phi += 2*M_PI;
      float theta = std::acos(dir.y);
      float u = phi/(2*M_PI);
      float v = theta/M_PI;
      return Vec2(u, v);
    };

    float terrainDist(const Vec3& p) const {
      Vec2 uv = approxUV(p);
      float sphereDist = (p - sphere->center).length() - sphere->radius;
      return sphereDist - heightMap(uv);
    };

  public:
    Terrain(const std::shared_ptr<Material>& material, const Vec3& origin, float radius) : Shape(material) {
      sphere = std::make_shared<Sphere>(material, origin, radius);
    };

    bool intersect(const Ray& ray, Hit& res) const {
      //二分探索
      if(sphere->intersect(ray, res)) {
        Vec3 startPoint = ray.origin;
        Vec3 endPoint = res.hitPos;
        for(int i = 0; i < 1000; i++) {
          Vec3 midPoint = (startPoint + endPoint)/2;
          float tdist = terrainDist(midPoint);

          if(std::abs(tdist) < 0.001) {
            res.t = (midPoint - ray.origin).length();
            res.hitPos = midPoint;
            Vec2 uv = approxUV(midPoint);
            res.u = uv.x;
            res.v = uv.y;
            res.hitNormal = heightMapNormal(midPoint);
            return true;
            break;
          }

          if(tdist > 0) {
            startPoint = midPoint;
          }
          else {
            endPoint = midPoint;
          }
        }
        return false;
      }
      //レイマーチング
      else {
        bool hit = false;
        float t = 0;
        Vec3 p = ray(t);
        res.iteration = 0;
        for(int i = 0; i < 1000; i++) {
          float dist = terrainDist(p);
          res.iteration++;
          if(std::abs(dist) < 0.001) {
            hit = true;
            res.t = t;
            res.hitPos = p;
            Vec2 uv = approxUV(p);
            res.u = uv.x;
            res.v = uv.y;
            res.hitNormal = heightMapNormal(p);
            break;
          }
          t += dist/2;
          p = ray(t);
        }
        return hit;
      }
    };
};
#endif
