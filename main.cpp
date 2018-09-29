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


const int samples = 1;
const int volume_samples = 10;
const int scattering_depth = 3;
const float R = 6360;
const float R_atmos = 6420;


const Vec3 sunDir = normalize(Vec3(0, -1, 0));
const RGB sunColor = RGB(100);


const float rayleigh_scaleheight = 8;
const float mie_scaleheight = 2.4;
const RGB beta_rayleigh = RGB(3.8e-3, 13.5e-3, 33.1e-3);
const RGB beta_mie = RGB(21e-3);


float rayleigh_beta_lambda(float height, float lambda) {
  return 8*std::pow(M_PI, 3.0)*std::pow(std::pow(1.000292, 2.0) - 1, 2.0)/(3*1.2*std::pow(lambda, 4.0)) * std::exp(-height/8.0);
}
RGB rayleigh_beta(float height) {
  return RGB(rayleigh_beta_lambda(height, 440*1e-9), rayleigh_beta_lambda(height, 550*1e-9), rayleigh_beta_lambda(height, 680*1e-9));
}


RGB Tr(const Vec3& p1, const Vec3& p2) {
  Ray lightRay(p1, normalize(p2 - p1));
  float rayleigh_optical_depth_light = 0;
  float mie_optical_depth_light = 0;
  float ds_light = (p2 - p1).length()/volume_samples;
  for(int j = 0; j < volume_samples; j++) {
    Vec3 p_light = lightRay(ds_light * (j + 1));
    float h_light = p_light.length() - R;
    rayleigh_optical_depth_light += std::exp(-h_light/rayleigh_scaleheight) * ds_light;
    mie_optical_depth_light += std::exp(-h_light/mie_scaleheight) * ds_light;
  }
  return exp(-(beta_rayleigh * rayleigh_optical_depth_light + beta_mie * mie_optical_depth_light));
}


float rayleigh_phase_function(const Vec3& wo, const Vec3& wi) {
  float mu = dot(wo, wi);
  return 3.0/(16.0*M_PI) * (1.0 + mu*mu);
}
float mie_phase_function(const Vec3& wo, const Vec3& wi, float g) {
  double mu = dot(wo, wi);
  return 3.0/(8.0*M_PI) * ((1.0 - g*g)*(1.0 + mu*mu))/((2.0 + g*g)*std::pow(1.0 + g*g - 2.0*g*mu, 1.5));
}


RGB Li(const Ray& _ray, const Scene& scene, Sampler& sampler) {
  Ray ray = _ray;
  RGB L;

  bool insideAtmos = ray.origin.length() - R_atmos < 0;
  if(!insideAtmos) {
    Hit res;
    if(scene.intersect(ray, res) && res.hitShape->type == "atmos") {
      ray = Ray(res.hitPos, ray.direction);
    }
    else {
      return RGB(0);
    }
  }

  Hit res;
  if(scene.intersect(ray, res)) {
    float ds = res.t/volume_samples;
    for(int i = 0; i < volume_samples; i++) {
      Vec3 p = ray(ds * (i + 1));
      RGB beta(1);
      RGB rayleigh_coeff(1);
      RGB mie_coeff(1);
      Vec3 wo = -ray.direction;
      for(int j = 0; j < scattering_depth; j++) {
        if(j != 0) {
          Vec3 p_prev = p;
          p = p_prev + ds*sampleSphere(sampler.getNext2D());
          beta *= 4*M_PI * Tr(p_prev, p);
          Vec3 wi = normalize(p - p_prev);
          rayleigh_coeff *= rayleigh_phase_function(wo, wi);
          mie_coeff *= mie_phase_function(wo, normalize(p - p_prev), 0.76);
          wo = -wi;
        }
        float h = p.length() - R;
        if(h < 0) break;
        float h_rayleigh = std::exp(-h/rayleigh_scaleheight) * ds;
        float h_mie = std::exp(-h/mie_scaleheight) * ds;
        RGB tr = Tr(ray.origin, p);
        beta *= tr;
        
        //Direct Light
        Ray lightRay = Ray(p, sunDir);
        Hit light_res;
        if(!scene.intersect(lightRay, light_res) || light_res.hitShape->type == "ground") break;
        rayleigh_coeff *= beta_rayleigh * h_rayleigh;
        mie_coeff *= beta_mie * h_mie;
        RGB tr_light = Tr(p, light_res.hitPos);
        L += (rayleigh_coeff * rayleigh_phase_function(wo, sunDir) + mie_coeff * mie_phase_function(wo, sunDir, 0.76)) * beta * tr_light * sunColor;
      }
    }

    if(res.hitShape->type == "ground") {
      auto hitMaterial = res.hitShape->material;
      float cos = std::max(dot(res.hitNormal, sunDir), 0.0);
      Ray lightRay = Ray(res.hitPos, sunDir);
      Hit light_res;
      scene.intersect(lightRay, light_res);

      if(light_res.hitShape->type == "atmos") {
        RGB tr= Tr(ray.origin, res.hitPos) * Tr(res.hitPos, light_res.hitPos);
        L += tr * hitMaterial->f(res, -ray.direction, sunDir) * cos * sunColor;
      }
    }
  }
  return L;
}


int main() {
  Film film(512, 512);
  Camera cam(Vec3(0, 0, -R - 1), normalize(Vec3(0, 1, 0)));

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

  for(int k = 0; k < samples; k++) {
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
