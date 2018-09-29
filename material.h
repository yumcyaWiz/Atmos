#ifndef MATERIAL_H
#define MATERIAL_H
#include "vec3.h"
#include "hit.h"
#include "sampler.h"
class Material {
  public:
    virtual RGB f(const Hit& res, const Vec3& wo, const Vec3& wi) const = 0;
    virtual RGB sample(const Hit& res, const Vec3& wo, Vec3& wi, Sampler& sampler, double& pdf) const = 0;
};
#endif
