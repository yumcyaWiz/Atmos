#ifndef VEC3_H
#define VEC3_H
#include <iostream>
#include <cmath>
class Vec3 {
  public:
    float x;
    float y;
    float z;

    Vec3() { x = y = z = 0; };
    Vec3(float _x) { x = y = z = _x; };
    Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {};

    Vec3 operator-() const {
      return Vec3(-x, -y, -z);
    };
    Vec3& operator+=(const Vec3& v) {
      x += v.x; y += v.y; z += v.z;
      return *this;
    };

    float length() const {
      return std::sqrt(x*x + y*y + z*z);
    };
    float length2() const {
      return x*x + y*y + z*z;
    };
};
inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}
inline Vec3 operator+(const Vec3& v1, float k) {
  return Vec3(v1.x + k, v1.y + k, v1.z + k);
}
inline Vec3 operator+(float k, const Vec3& v2) {
  return Vec3(k + v2.x, k + v2.y, k + v2.z);
}

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}
inline Vec3 operator-(const Vec3& v1, float k) {
  return Vec3(v1.x - k, v1.y - k, v1.z - k);
}
inline Vec3 operator-(float k, const Vec3& v2) {
  return Vec3(k - v2.x, k - v2.y, k - v2.z);
  }

inline Vec3 operator*(const Vec3 v1, const Vec3& v2) {
  return Vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline Vec3 operator*(const Vec3& v1, float k) {
  return Vec3(v1.x * k, v1.y * k, v1.z * k);
}
inline Vec3 operator*(float k, const Vec3& v2) {
  return Vec3(k * v2.x, k * v2.y, k * v2.z);
}

inline Vec3 operator/(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}
inline Vec3 operator/(const Vec3& v1, float k) {
  return Vec3(v1.x / k, v1.y / k, v1.z / k);
}
inline Vec3 operator/(float k, const Vec3& v2) {
  return Vec3(k / v2.x, k / v2.y, k / v2.z);
}


inline std::ostream& operator<<(std::ostream& stream, const Vec3& v) {
  stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return stream;
}


inline float dot(const Vec3& v1, const Vec3& v2) {
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}


inline Vec3 normalize(const Vec3& v) {
  return v/v.length();
}


inline Vec3 pow(const Vec3& v, float n) {
  return Vec3(std::pow(v.x, n), std::pow(v.y, n), std::pow(v.z, n));
}


inline void orthonormalBasis(const Vec3& n, Vec3& vx, Vec3& vz) {
  if(std::abs(n.x) > 0.9) vx = Vec3(0, 1, 0);
  else vx = Vec3(1, 0, 0);

  vx = normalize(vx - n*dot(vx, n));
  vz = cross(n, vx);
}


using RGB = Vec3;
#endif
