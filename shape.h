#ifndef SHAPE_H
#define SHAPE_H
#include <memory>
#include "material.h"
class Shape {
  public:
    std::shared_ptr<Material> material;
    std::string type;

    Shape(const std::shared_ptr<Material>& _material, const std::string& _type) : material(_material), type(_type) {};

    virtual bool intersect(const Ray& ray, Hit& res) const = 0;
};
#endif
