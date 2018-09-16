#ifndef FILM_H
#define FILM_H
#include <iostream>
#include <fstream>
#include <string>
#include "vec3.h"
#include "util.h"
class Film {
  private:
    struct Pixel {
      RGB color;
      int samples;
      Pixel() : color(RGB(0)), samples(0) {};
    };

  public:
    int width;
    int height;
    Pixel* data;

    Film(int _width, int _height) : width(_width), height(_height) {
      data = new Pixel[width*height];
    };
    ~Film() {
      delete[] data;
    };

    void addSample(int i, int j, const RGB& L) {
      if(i < 0 || i >= width || j < 0 || j >= height) {
        std::cerr << "Invalid indexes" << std::endl;
        std::exit(1);
      }
      data[i + width*j].color += L;
      data[i + width*j].samples++;
    };

    void ppm_output(const std::string& filename) const {
      std::ofstream file(filename);
      file << "P3" << std::endl;
      file << width << " " << height << std::endl;
      file << "255" << std::endl;
      for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
          RGB col = pow(data[i + width*j].color/data[i + width*j].samples, 1/2.2);
          int r = clamp(int(255*col.x), 0, 255);
          int g = clamp(int(255*col.y), 0, 255);
          int b = clamp(int(255*col.z), 0, 255);
          file << r << " " << g << " " << b << std::endl;
        }
      }
      file.close();
    };
};
#endif
