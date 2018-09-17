#ifndef ACCEL_H
#define ACCEL_H
#include <vector>
#include <memory>
#include "ray.h"
#include "hit.h"
template <typename T>
class Accel {
  public:
    std::vector<std::shared_ptr<T>> prims;

    Accel() {};
    Accel(const std::vector<std::shared_ptr<T>>& _prims) : prims(_prims) {};

    void add(const std::shared_ptr<T>& p) {
      prims.push_back(p);
    };

    bool intersect(const Ray& ray, Hit& res) const {
      bool hit = false;
      for(auto p : prims) {
        Hit res_temp;
        if(p->intersect(ray, res_temp)) {
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
