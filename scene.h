#ifndef SCENE_H
#define SCENE_H
#include <memory>
#include "ray.h"
#include "hit.h"
#include "shape.h"
#include "accel.h"
class Scene {
  public:
    std::shared_ptr<Accel<Shape>> accel;

    Scene(const std::vector<std::shared_ptr<Shape>>& shapes) {
      accel = std::make_shared<Accel<Shape>>(shapes);
    };

    bool intersect(const Ray& ray, Hit& res) const {
      return accel->intersect(ray, res);
    };
};
#endif
