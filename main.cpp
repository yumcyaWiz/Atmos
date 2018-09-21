#include "vec3.h"
#include "ray.h"
#include "shapes/sphere.h"
#include "shapes/terrain.h"
#include "film.h"
#include "camera.h"
#include "scene.h"
#include "textures/imageTexture.h"
#include "materials/lambert.h"


const float R = 10;
const float R_atmos = 6381 * 1000;


const Vec3 sunDir = normalize(Vec3(1, 1, -1));


RGB Li(const Ray& ray, const Scene& scene) {
  Hit res;
  if(scene.intersect(ray, res)) {
    return M_PI*res.hitShape->material->f(res, -ray.direction, sunDir) * std::max(dot(res.hitNormal, sunDir), 0.0f);
  }
  else {
    return RGB(0);
  }
}


int main() {
  Film film(512, 512);
  Camera cam(Vec3(0, 0, -2*R), normalize(Vec3(0, 0, 1)));

  auto tex = std::make_shared<ImageTexture>("earth.jpg");
  auto mat = std::make_shared<Lambert>(tex);

  std::vector<std::shared_ptr<Shape>> shapes;
  auto earth = std::make_shared<Sphere>(mat, Vec3(0, 0, 0), R);
  auto atmos = std::make_shared<Sphere>(mat, Vec3(0, 0, 0), R_atmos);
  auto terrain = std::make_shared<Terrain>(mat, Vec3(0, 0, 0), R);
  shapes.push_back(earth);

  Scene scene(shapes);

#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < film.width; i++) {
    for(int j = 0; j < film.height; j++) {
      float u = (2.0*i - film.width)/film.width;
      float v = (2.0*j - film.height)/film.width;
      Ray ray = cam.getRay(u, v);
      film.addSample(i, j, Li(ray, scene));
    }
  }
  film.ppm_output("output.ppm");

  return 0;
}
