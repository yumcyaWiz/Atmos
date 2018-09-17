#ifndef SHAPE_H
#define SHAPE_H
class Shape {
  public:
    virtual bool intersect(const Ray& ray, Hit& res) const = 0;
};
#endif
