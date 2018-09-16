#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "film.h"
#include "camera.h"


int main() {
  Film film(512, 512);
  Camera cam(Vec3(0, 0, 0), Vec3(0, 0, 1));

  Sphere sphere(Vec3(0, 0, 3), 1.0);

#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < film.width; i++) {
    for(int j = 0; j < film.height; j++) {
      float u = (2.0*i - film.width)/film.width;
      float v = (2.0*j - film.height)/film.width;
      Ray ray = cam.getRay(u, v);
      
      Hit res;
      if(sphere.intersect(ray, res)) {
        film.addSample(i, j, (res.hitNormal + 1)/2);
      }
    }
  }
  film.ppm_output("output.ppm");

  return 0;
}
