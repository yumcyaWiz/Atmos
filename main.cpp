#include "vec3.h"
#include "ray.h"
#include "shapes/sphere.h"
#include "shapes/terrain.h"
#include "film.h"
#include "camera.h"
#include "scene.h"
#include "textures/imageTexture.h"
#include "materials/lambert.h"
#include "samplers/mt.h"


const float R = 6371;
const float R_atmos = 6381;


const Vec3 sunDir = normalize(Vec3(1, 1, -1));
const RGB sunColor = RGB(1, 1, 1);


const float rayleigh_scaleheight = 7.4;
const float mie_scaleheight = 1.2;
const RGB beta_rayleigh = RGB(3.3e-6, 13.6e-6, 33.1e-6);
const RGB beta_mie = RGB(21e-6);
float rayleigh_phase_function(const Vec3& wo, const Vec3& wi) {
  float mu = dot(wo, wi);
  return 3.0/(16.0*M_PI) * (1.0 + mu*mu);
}
float mie_phase_function(const Vec3& wo, const Vec3& wi, float g) {
  double mu = dot(wo, wi);
  return 3.0/(8.0*M_PI) * ((1.0 - g*g)*(1.0 + mu*mu))/((2.0 + g*g)*std::pow(1.0 + g*g - 2.0*g*mu, 1.5));
}


const int maxDepth = 10000;
const float ds = 1;
RGB Li(const Ray& _ray, const Scene& scene, Sampler& sampler) {
  RGB L;
  RGB T(1);

  Ray ray = _ray;
  for(int depth = 0; depth < maxDepth; depth++) {
    if(isZero(T)) break;

    Hit res;
    if(scene.intersect(ray, res)) {
      float h = ray.origin.length() - R;
      Vec3 wo = -ray.direction;
      float cos = std::max(dot(wo, sunDir), 0.0f);

      //地面に当たった場合
      if(res.t < ds && res.hitShape->type == "ground") {
        auto hitMaterial = res.hitShape->material;
        L += T * M_PI*hitMaterial->f(res, wo, sunDir)*cos*sunColor;
        break;
      }
      //大気圏外から大気に当たった場合
      else if(h > (R_atmos - R) && res.hitShape->type == "atmos") {
        //Next Ray
        ray = Ray(res.hitPos, ray.direction);
      }
      else {
        Vec3 p = ray(ds);
        h = p.length() - R;

        float optical_depth_rayleigh = std::exp(-h/rayleigh_scaleheight) * ds;
        RGB tau = beta_rayleigh * optical_depth_rayleigh;
        T *= exp(-tau);

        //Direct Light Sampling
        Ray shadowRay(p, sunDir);
        Hit shadow_res;
        bool hit = scene.intersect(shadowRay, shadow_res);
        if(!hit || shadow_res.hitShape->type == "atmos") {
          L += T * rayleigh_phase_function(wo, shadowRay.direction) * sunColor;
        }

        //NextRay
        Vec3 wi = sampleSphere(sampler.getNext2D());
        ray = Ray(p, wi);
        T *= rayleigh_phase_function(wo, wi);
      }
    }
    else {
      break;
    }
  }
  return L;
}


int main() {
  Film film(512, 512);
  Camera cam(Vec3(R + 1, 0, 0), normalize(Vec3(0, 0, 1)));

  auto tex = std::make_shared<ImageTexture>("earth2.jpg");
  auto mat = std::make_shared<Lambert>(tex);

  std::vector<std::shared_ptr<Shape>> shapes;
  auto earth = std::make_shared<Sphere>(mat, "ground", Vec3(0, 0, 0), R);
  auto atmos = std::make_shared<Sphere>(mat, "atmos", Vec3(0, 0, 0), R_atmos);
  auto terrain = std::make_shared<Terrain>(mat, "ground", Vec3(0, 0, 0), R);
  shapes.push_back(earth);
  shapes.push_back(atmos);

  Scene scene(shapes);
  Mt mt;

  for(int k = 0; k < 1; k++) {
#pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < film.width; i++) {
      for(int j = 0; j < film.height; j++) {
        float u = (2.0*i + mt.getNext() - film.width)/film.width;
        float v = (2.0*j + mt.getNext() - film.height)/film.width;
        Ray ray = cam.getRay(u, v);
        film.addSample(i, j, Li(ray, scene, mt));
      }
    }
  }
  film.ppm_output("output.ppm");

  return 0;
}
