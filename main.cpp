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
const int volume_samples = 10;
const int scattering_depth = 3;
const double R = 6360 * 1000;
const double R_atmos = 6420 * 1000;


const Vec3 sunDir = normalize(Vec3(0, 1, -0.1));
const RGB sunColor = RGB(20);


const double rayleigh_scaleheight = 8000;
const double mie_scaleheight = 1200;
const RGB beta_rayleigh = RGB(3.8e-6, 13.5e-6, 33.1e-6);
const RGB beta_mie = RGB(21e-6);


double rayleigh_beta_lambda(double height, double lambda) {
  return 8*std::pow(M_PI, 3.0)*std::pow(std::pow(1.000292, 2.0) - 1, 2.0)/(3*1.2*std::pow(lambda, 4.0)) * std::exp(-height/8.0);
}
RGB rayleigh_beta(double height) {
  return RGB(rayleigh_beta_lambda(height, 440*1e-9), rayleigh_beta_lambda(height, 550*1e-9), rayleigh_beta_lambda(height, 680*1e-9));
}


RGB Tr(const Vec3& p1, const Vec3& p2) {
  Ray lightRay(p1, normalize(p2 - p1));
  double rayleigh_optical_depth = 0;
  double mie_optical_depth = 0;
  double ds_light = (p2 - p1).length()/volume_samples;
  for(int j = 0; j < volume_samples; j++) {
    Vec3 p_light = lightRay(ds_light * (j + 1));
    double h_light = p_light.length() - R;
    rayleigh_optical_depth += std::exp(-h_light/rayleigh_scaleheight) * ds_light;
    mie_optical_depth += std::exp(-h_light/mie_scaleheight) * ds_light;
  }
  return exp(-(beta_rayleigh * rayleigh_optical_depth + 1.1 * beta_mie * mie_optical_depth));
}


double rayleigh_phase_function(const Vec3& wo, const Vec3& wi) {
  double mu = dot(wo, -wi);
  return 3.0/(16.0*M_PI) * (1.0 + mu*mu);
}
double mie_phase_function(const Vec3& wo, const Vec3& wi, double g) {
  double mu = dot(wo, -wi);
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
    double ds = res.t/volume_samples;
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
          mie_coeff *= mie_phase_function(wo, wi, 0.76);
          wo = -wi;
        }
        else {
          beta *= Tr(ray.origin, p);
        }
        double h = p.length() - R;
        if(h < 0) break;
        double h_rayleigh = std::exp(-h/rayleigh_scaleheight) * ds;
        double h_mie = std::exp(-h/mie_scaleheight) * ds;
        
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
      double cos = std::max(dot(res.hitNormal, sunDir), 0.0);
      Ray lightRay = Ray(res.hitPos + res.hitNormal, sunDir);
      /*
      Hit light_res;
      scene.intersect(lightRay, light_res);

      if(light_res.hitShape->type == "atmos") {
        RGB tr = Tr(ray.origin, res.hitPos) * Tr(res.hitPos, light_res.hitPos);
        L += tr * hitMaterial->f(res, -ray.direction, sunDir) * cos * sunColor;
      }
      */
      L += hitMaterial->f(res, -ray.direction, sunDir) * cos * Tr(ray.origin, res.hitPos) * Li(lightRay, scene, sampler);
    }
  }
  return L;
}


int main() {
  Film film(2138, 1536);
  Camera cam(Vec3(0, 0, -R - 0.1*1000), normalize(Vec3(0, 1, 0)));

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
        double u = (2.0*i + mt.getNext() - film.width)/film.width;
        double v = (2.0*j + mt.getNext() - film.height)/film.width;
        Ray ray = cam.getRay(u, v);
        film.addSample(i, j, Li(ray, scene, mt));
      }
    }
    std::cout << double(k)/samples * 100 << "%" << std::endl;
  }
  film.ppm_output("output.ppm");

  return 0;
}
