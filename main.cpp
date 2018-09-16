#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "film.h"


int main() {
  Film film(512, 512);

#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < film.width; i++) {
    for(int j = 0; j < film.height; j++) {
      float u = (2.0*i - film.width)/film.width;
      float v = (2.0*j - film.height)/film.width;
      film.addSample(i, j, RGB((u + 1)/2, (v + 1)/2, 1));
    }
  }
  film.ppm_output("output.ppm");

  return 0;
}
