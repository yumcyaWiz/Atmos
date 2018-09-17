#include "vec3.h"
#include "ray.h"
#include "shapes/sphere.h"
#include "shapes/terrain.h"
#include "film.h"
#include "camera.h"
#include "scene.h"


const float R = 1000;
const float R_atmos = 6381 * 1000;


int main() {
  Film film(512, 512);
  Camera cam(Vec3(0, R + 10, 100), normalize(Vec3(0, 0, 1)));

  std::vector<std::shared_ptr<Shape>> shapes;
  auto earth = std::make_shared<Sphere>(Vec3(0, 0, 0), R);
  auto atmos = std::make_shared<Sphere>(Vec3(0, 0, 0), R_atmos);
  auto terrain = std::make_shared<Terrain>(Vec3(0, 0, 0), R);
  shapes.push_back(terrain);

  Scene scene(shapes);

#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < film.width; i++) {
    for(int j = 0; j < film.height; j++) {
      float u = (2.0*i - film.width)/film.width;
      float v = (2.0*j - film.height)/film.width;
      Ray ray = cam.getRay(u, v);
      Hit res;
      if(scene.intersect(ray, res)) {
        film.addSample(i, j, RGB(dot(res.hitNormal, Vec3(0, 1, 0))));
      }
    }
  }
  film.ppm_output("output.ppm");

  return 0;
}
