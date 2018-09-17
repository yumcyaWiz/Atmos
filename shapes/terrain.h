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
      return std::sin(100*u.x)*std::sin(100*u.y) + 1;
    };
    Vec3 heightMapNormal(const Vec3& p) const {
      Vec3 n = normalize(p - sphere->center);
      Vec2 u = approxUV(p);
      float du = 0.1*(heightMap(u + Vec2(0.001, 0)) - heightMap(u - Vec2(0.001, 0)))/0.002;
      float dv = 0.1*(heightMap(u + Vec2(0, 0.001)) - heightMap(u - Vec2(0, 0.001)))/0.002;
      return normalize(Vec3(du, 1, dv));
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
    Terrain(const Vec3& origin, float radius) {
      sphere = std::make_shared<Sphere>(origin, radius);
    };

    bool intersect(const Ray& ray, Hit& res) const {
      //二分探索
      if(sphere->intersect(ray, res)) {
        Vec3 startPoint = ray.origin;
        Vec3 endPoint = res.hitPos;
        for(int i = 0; i < 1000; i++) {
          Vec3 midPoint = (startPoint + endPoint)/2;
          float tdist = terrainDist(midPoint);

          if(std::abs(tdist) < 0.1) {
            res.t = (midPoint - ray.origin).length();
            res.hitPos = midPoint;
            Vec2 uv = approxUV(res.hitPos);
            res.u = uv.x;
            res.v = uv.y;
            res.hitNormal = heightMapNormal(res.hitPos);
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
          if(std::abs(dist) < 0.1) {
            hit = true;
            res.t = t;
            res.hitPos = p;
            Vec2 uv = approxUV(p);
            res.u = uv.x;
            res.v = uv.y;
            res.hitNormal = heightMapNormal(p);
            break;
          }
          t += dist;
          p = ray(t);
        }
        return hit;
      }
    };
};
#endif
