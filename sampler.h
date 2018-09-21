#ifndef SAMPLER_H
#define SAMPLER_H
#include "vec2.h"
class Sampler {
  public:
    virtual float getNext() const = 0;
    virtual Vec2 getNext2D() const = 0;
};
#endif
