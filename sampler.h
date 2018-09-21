#ifndef SAMPLER_H
#define SAMPLER_H
#include "vec2.h"
class Sampler {
  public:
    virtual float getNext() = 0;
    virtual Vec2 getNext2D() = 0;
};
#endif
