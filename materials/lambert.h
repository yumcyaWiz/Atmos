#ifndef LAMBERT_H
#define LAMBERT_H
#include "../material.h"
class Lambert : public Material {
  public:
    std::shared_ptr<Texture> reflectance;
    
    Lambert(const std::shared_ptr<Texture>& _reflectance) : reflectance(_reflectance) {};

    RGB f(const Hit& res, const Vec3& wo, const Vec3& wi) const {
      return reflectance->getColor(res)/M_PI;
    };
    RGB sample(const Hit& res, const Vec3& wo, Vec3& wi, Sampler& sampler, float& pdf) const {
      return reflectance->getColor(res)/M_PI;
    };
};
#endif
