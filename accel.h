#ifndef ACCEL_H
#define ACCEL_H
#include <vector>
#include <memory>
#include "ray.h"
#include "hit.h"
#include "sphere.h"
class Accel {
  public:
    std::vector<std::shared_ptr<Sphere>> spheres;

    Accel() {};
    Accel(const std::vector<std::shared_ptr<Sphere>>& _spheres) : spheres(_spheres) {};

    void add(const std::shared_ptr<Sphere>& s) {
      spheres.push_back(s);
    };

    bool intersect(const Ray& ray, Hit& res) const {
      bool hit = false;
      for(auto s : spheres) {
        Hit res_temp;
        if(s->intersect(ray, res_temp)) {
          if(res_temp.t < res.t) {
            hit = true;
            res = res_temp;
          }
        }
      }
      return hit;
    };
};
#endif
