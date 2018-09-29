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


const int samples = 100;
const float R = 6360;
const float R_atmos = 6420;


const Vec3 sunDir = normalize(Vec3(0, 0, -1));
const RGB sunColor = RGB(5);


const float rayleigh_scaleheight = 8.0;
const float mie_scaleheight = 1.2;
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
  float ds_light = (p2 - p1).length()/samples;
  for(int j = 0; j < samples; j++) {
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
    float ds = res.t/samples;
    for(int i = 0; i < samples; i++) {
      Vec3 p = ray(ds * (i + 1));
      float h = p.length() - R;
      float h_rayleigh = std::exp(-h/rayleigh_scaleheight) * ds;
      float h_mie = std::exp(-h/mie_scaleheight) * ds;
      RGB tr = Tr(ray.origin, p);
      
      //一次散乱
      Ray lightRay = Ray(p, sunDir);
      Hit light_res;
      if(!scene.intersect(lightRay, light_res) || light_res.hitShape->type == "ground") break;
      RGB rayleigh_coeff = beta_rayleigh * h_rayleigh;
      RGB mie_coeff = beta_mie * h_mie;
      RGB tr_light = Tr(p, light_res.hitPos);
      L += (rayleigh_coeff * rayleigh_phase_function(-ray.direction, sunDir) + mie_coeff * mie_phase_function(-ray.direction, sunDir, 0.76)) * tr * tr_light * sunColor;

      //二次散乱
      Vec3 p2 = p + ds*sampleSphere(sampler.getNext2D());
      lightRay = Ray(p2, sunDir);
      light_res = Hit();
      if(!scene.intersect(lightRay, light_res) || light_res.hitShape->type == "ground") {
        h = p2.length() - R;
        float rh = std::exp(-h/rayleigh_scaleheight);
        float mh = std::exp(-h/mie_scaleheight);
        RGB rayleigh_coeff2 = beta_rayleigh * rh;
        RGB mie_coeff2 = beta_mie * mh;
        RGB tr2 = Tr(p, p2) * Tr(p2, light_res.hitPos);
        L += (rayleigh_coeff*rayleigh_coeff2 * rayleigh_phase_function(-ray.direction, sunDir)*rayleigh_phase_function(normalize(p - p2), sunDir) + mie_coeff*mie_coeff2 * mie_phase_function(-ray.direction, sunDir, 0.76)*mie_phase_function(normalize(p - p2), sunDir, 0.76)) * tr * tr2 * sunColor;
      }
    }

    if(res.hitShape->type == "ground") {
      auto hitMaterial = res.hitShape->material;
      float cos = std::max(dot(res.hitNormal, sunDir), 0.0f);
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
  Camera cam(Vec3(0, 0, -R - 10), normalize(Vec3(0, 1, 0)));

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
