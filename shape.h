#ifndef SHAPE_H
#define SHAPE_H
#include <memory>
#include "material.h"
class Shape {
  public:
    std::shared_ptr<Material> material;

    Shape(const std::shared_ptr<Material>& _material) : material(_material) {};

    virtual bool intersect(const Ray& ray, Hit& res) const = 0;
};
#endif
