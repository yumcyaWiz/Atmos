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


const float R = 6360;
const float R_atmos = 6420;


const Vec3 sunDir = normalize(Vec3(0, 0, -1));
const RGB sunColor = RGB(5);


const float rayleigh_scaleheight = 8.0;
const float mie_scaleheight = 1.2;
const RGB beta_rayleigh = RGB(3.8e-3, 13.5e-3, 33.1e-3);
const RGB beta_mie = RGB(21e-3);
float rayleigh_phase_function(const Vec3& wo, const Vec3& wi) {
  float mu = dot(wo, wi);
  return 3.0/(16.0*M_PI) * (1.0 + mu*mu);
}
float mie_phase_function(const Vec3& wo, const Vec3& wi, float g) {
  double mu = dot(wo, wi);
  return 3.0/(8.0*M_PI) * ((1.0 - g*g)*(1.0 + mu*mu))/((2.0 + g*g)*std::pow(1.0 + g*g - 2.0*g*mu, 1.5));
}


const int samples = 100;
RGB Li(const Ray& _ray, const Scene& scene, Sampler& sampler) {
  Ray ray = _ray;
  RGB L;
  Hit res;
  bool insideAtmos = ray.origin.length() - R_atmos < 0;
  if(!insideAtmos) {
    if(scene.intersect(ray, res) && res.hitShape->type == "atmos") {
      ray = Ray(res.hitPos, ray.direction);
    }
    else {
      return RGB(0);
    }
  }

  if(scene.intersect(ray, res)) {
    float ds = res.t/samples;
    float rayleigh_optical_depth = 0;
    for(int i = 0; i < samples; i++) {
      Vec3 p = ray(ds * (i + 1));
      float h = p.length() - R;
      float h_rayleigh = std::exp(-h/rayleigh_scaleheight) * ds;
      rayleigh_optical_depth += h_rayleigh;

      Ray lightRay = Ray(p, sunDir);
      Hit light_res;
      if(!scene.intersect(lightRay, light_res)) {
        break;
      }
      
      float rayleigh_optical_depth_light = 0;
      float ds_light = light_res.t/samples;
      for(int j = 0; j < samples; j++) {
        Vec3 p_light = lightRay(ds_light * (j + 1));
        float h_light = p_light.length() - R;
        rayleigh_optical_depth_light += std::exp(-h_light/rayleigh_scaleheight) * ds_light;
      }

      RGB trans = exp(-beta_rayleigh * (rayleigh_optical_depth + rayleigh_optical_depth_light));
      RGB coeff_s = beta_rayleigh * h_rayleigh;
      L += coeff_s * rayleigh_phase_function(-ray.direction, sunDir) * trans * sunColor;
    }

    if(res.hitShape->type == "ground") {
      auto hitMaterial = res.hitShape->material;
      float cos = std::max(dot(-ray.direction, sunDir), 0.0f);
      Ray lightRay = Ray(res.hitPos, sunDir);
      Hit light_res;
      scene.intersect(lightRay, light_res);
      float rayleigh_optical_depth_light = 0;
      float ds_light = light_res.t/samples;
      for(int j = 0; j < samples; j++) {
        Vec3 p_light = lightRay(ds_light * (j + 1));
        float h_light = p_light.length() - R;
        rayleigh_optical_depth_light += std::exp(-h_light/rayleigh_scaleheight) * ds_light;
      }
      RGB trans = exp(-beta_rayleigh * (rayleigh_optical_depth + rayleigh_optical_depth_light));
      L += trans * hitMaterial->f(res, -ray.direction, sunDir) * cos * sunColor;
    }
  }
  return L;
}


int main() {
  Film film(512, 512);
  Camera cam(Vec3(0, 0, -2*R), normalize(Vec3(0, 0, 1)));

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
