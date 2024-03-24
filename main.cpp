#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <strstream>
#include <utility>
#include <vector>

#define RESERVE_VERTS 64
#define RESERVE_TRIS RESERVE_VERTS * 2
#define RESERVE_RAST_TRIS RESERVE_TRIS / 4
#define NEAR_PLANE 0.1f

struct vec2d {
  float u;
  float v;
  float w;
  vec2d() : u(0.0f), v(0.0f), w(1.0f) {}
  vec2d(const float x) : u(x), v(x), w(1.0f) {}
  vec2d(const float u, const float v) : u(u), v(v), w(1.0f) {}

  vec2d &operator/=(const float k) {
    this->u /= k;
    this->v /= k;
    return *this;
  }
};

struct vec3d {
  float x, y, z;
  float w;

  vec3d() : x(0.0f), y(0.0f), z(0.0f), w(1.0f) {}

  vec3d(const float x) : x(x), y(x), z(x), w(1.0f) {}

  vec3d(const float x, const float y, const float z)
      : x(x), y(y), z(z), w(1.0f) {}

  vec3d(const float x, const float y, const float z, const float w)
      : x(x), y(y), z(z), w(w) {}

  vec3d &operator+=(const vec3d &other) {
    this->x += other.x;
    this->y += other.y;
    this->z += other.z;
    return *this;
  }

  vec3d &operator-=(const vec3d &other) {
    this->x -= other.x;
    this->y -= other.y;
    this->z -= other.z;
    return *this;
  }

  vec3d &operator*=(const vec3d &other) {
    this->x *= other.x;
    this->y *= other.y;
    this->z *= other.z;
    return *this;
  }

  vec3d &operator*=(float k) {
    this->x *= k;
    this->y *= k;
    this->z *= k;
    return *this;
  }

  vec3d &operator/=(const vec3d &other) {
    this->x /= other.x;
    this->y /= other.y;
    this->z /= other.z;
    return *this;
  }

  vec3d &operator/=(float k) {
    this->x /= k;
    this->y /= k;
    this->z /= k;
    return *this;
  }

  float dot_product(const vec3d &other) const {
    return this->x * other.x + this->y * other.y + this->z * other.z;
  }

  float length() const { return sqrtf(this->dot_product(*this)); }

  void normalize() { *this /= this->length(); }

  vec3d cross_product(const vec3d &other) const {
    return {
        (this->y * other.z - this->z * other.y),
        (this->z * other.x - this->x * other.z),
        (this->x * other.y - this->y * other.x),
    };
  }
};

vec3d operator-(vec3d lhs, const vec3d &rhs) {
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  lhs.z -= rhs.z;
  return lhs;
}

vec3d operator+(vec3d lhs, const vec3d &rhs) {
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  lhs.z += rhs.z;
  return lhs;
}

vec3d operator*(vec3d lhs, const float rhs) {
  lhs.x *= rhs;
  lhs.y *= rhs;
  lhs.z *= rhs;
  return lhs;
}

vec3d vector_intersect_plane(const vec3d &plane_p, vec3d &plane_n,
                             const vec3d &lineStart, const vec3d &lineEnd,
                             float &t) {
  plane_n.normalize();
  float plane_d = -plane_n.dot_product(plane_p);
  float ad = lineStart.dot_product(plane_n);
  float bd = lineEnd.dot_product(plane_n);
  t = (-plane_d - ad) / (bd - ad);
  vec3d lineStartToEnd = lineEnd - lineStart;
  vec3d lineToIntersect = lineStartToEnd * t;
  return lineStart + lineToIntersect;
}

struct mat4 {
  float m[4][4] = {{0.0f}};

  void identity() {
    this->m[0][0] = 1.0f;
    this->m[1][1] = 1.0f;
    this->m[2][2] = 1.0f;
    this->m[3][3] = 1.0f;
  }

  void rotation_x(float fAngleRad) {
    this->m[0][0] = 1.0f;
    this->m[3][3] = 1.0f;
    this->m[1][1] = cosf(fAngleRad);
    this->m[2][2] = cosf(fAngleRad);
    this->m[1][2] = sinf(fAngleRad);
    this->m[2][1] = -sinf(fAngleRad);
  }

  void rotation_y(float fAngleRad) {
    this->m[1][1] = 1.0f;
    this->m[3][3] = 1.0f;
    this->m[0][0] = cosf(fAngleRad);
    this->m[2][2] = cosf(fAngleRad);
    this->m[0][2] = sinf(fAngleRad);
    this->m[2][0] = -sinf(fAngleRad);
  }

  void rotation_z(float fAngleRad) {
    this->m[2][2] = 1.0f;
    this->m[3][3] = 1.0f;
    this->m[0][0] = cosf(fAngleRad);
    this->m[1][1] = cosf(fAngleRad);
    this->m[0][1] = sinf(fAngleRad);
    this->m[1][0] = -sinf(fAngleRad);
  }

  void translation(float x, float y, float z) {
    this->m[0][0] = 1.0f;
    this->m[1][1] = 1.0f;
    this->m[2][2] = 1.0f;
    this->m[3][3] = 1.0f;
    this->m[3][0] = x;
    this->m[3][1] = y;
    this->m[3][2] = z;
  }

  void projection(float fFovDegrees, float fAspectRatio, float fNear,
                  float fFar) {
    float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
    this->m[0][0] = fAspectRatio * fFovRad;
    this->m[1][1] = fFovRad;
    this->m[2][2] = fFar / (fFar - fNear);
    this->m[3][2] = (-fFar * fNear) / (fFar - fNear);
    this->m[2][3] = 1.0f;
    this->m[3][3] = 0.0f;
  }

  void point_at(const vec3d &pos, vec3d &target, vec3d &up) {
    target = target - pos;
    target.normalize();

    vec3d a = target * up.dot_product(target);
    up = up - a;
    up.normalize();

    vec3d newRight = up.cross_product(target);

    this->m[0][0] = newRight.x;
    this->m[0][1] = newRight.y;
    this->m[0][2] = newRight.z;
    this->m[0][3] = 0.0f;
    this->m[1][0] = up.x;
    this->m[1][1] = up.y;
    this->m[1][2] = up.z;
    this->m[1][3] = 0.0f;
    this->m[2][0] = target.x;
    this->m[2][1] = target.y;
    this->m[2][2] = target.z;
    this->m[2][3] = 0.0f;
    this->m[3][0] = pos.x;
    this->m[3][1] = pos.y;
    this->m[3][2] = pos.z;
    this->m[3][3] = 1.0f;
  }

  void quick_inverse() {
    mat4 matrix;
    matrix.m[0][0] = this->m[0][0];
    matrix.m[0][1] = this->m[1][0];
    matrix.m[0][2] = this->m[2][0];
    matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = this->m[0][1];
    matrix.m[1][1] = this->m[1][1];
    matrix.m[1][2] = this->m[2][1];
    matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = this->m[0][2];
    matrix.m[2][1] = this->m[1][2];
    matrix.m[2][2] = this->m[2][2];
    matrix.m[2][3] = 0.0f;
    matrix.m[3][0] =
        -(this->m[3][0] * matrix.m[0][0] + this->m[3][1] * matrix.m[1][0] +
          this->m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] =
        -(this->m[3][0] * matrix.m[0][1] + this->m[3][1] * matrix.m[1][1] +
          this->m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] =
        -(this->m[3][0] * matrix.m[0][2] + this->m[3][1] * matrix.m[1][2] +
          this->m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    *this = matrix;
  }

  mat4 operator*=(const mat4 &other) {
    mat4 matrix;
    for (int c = 0; c < 4; c++)
      for (int r = 0; r < 4; r++)
        matrix.m[r][c] =
            this->m[r][0] * other.m[0][c] + this->m[r][1] * other.m[1][c] +
            this->m[r][2] * other.m[2][c] + this->m[r][3] * other.m[3][c];
    return matrix;
  }
};

mat4 operator*(const mat4 &lhs, const mat4 &rhs) {
  mat4 matrix;
  for (int c = 0; c < 4; c++)
    for (int r = 0; r < 4; r++)
      matrix.m[r][c] = lhs.m[r][0] * rhs.m[0][c] + lhs.m[r][1] * rhs.m[1][c] +
                       lhs.m[r][2] * rhs.m[2][c] + lhs.m[r][3] * rhs.m[3][c];
  return matrix;
}

vec3d operator*(const vec3d &lhs, const mat4 &rhs) {
  vec3d v;
  v.x = lhs.x * rhs.m[0][0] + lhs.y * rhs.m[1][0] + lhs.z * rhs.m[2][0] +
        lhs.w * rhs.m[3][0];
  v.y = lhs.x * rhs.m[0][1] + lhs.y * rhs.m[1][1] + lhs.z * rhs.m[2][1] +
        lhs.w * rhs.m[3][1];
  v.z = lhs.x * rhs.m[0][2] + lhs.y * rhs.m[1][2] + lhs.z * rhs.m[2][2] +
        lhs.w * rhs.m[3][2];
  v.w = lhs.x * rhs.m[0][3] + lhs.y * rhs.m[1][3] + lhs.z * rhs.m[2][3] +
        lhs.w * rhs.m[3][3];
  return v;
}

struct triangle {
  vec3d p[3];
  vec2d t[3];

  olc::Pixel col;

  triangle() : p(), col(olc::BLACK) {}

  triangle(const vec3d &x, const vec3d &y, const vec3d &z) : p{x, y, z} {}

  triangle(const vec3d &x, const vec3d &y, const vec3d &z, const vec2d a,
           const vec2d b, const vec2d c)
      : p{x, y, z}, t{a, b, c} {}

  triangle(const vec3d &x, const vec3d &y, const vec3d &z,
           const olc::Pixel color)
      : p{x, y, z}, col(color) {}

#ifdef OPTIM_DEBUG
  triangle &operator=(const triangle &other) {
    p[0] = other.p[0];
    p[1] = other.p[1];
    p[2] = other.p[2];
    sym = other.sym;
    col = other.col;
    std::cout << "copied\n";
    return *this;
  }
#endif

  triangle &operator+=(const triangle &other) {
    this->p[0] += other.p[0];
    this->p[1] += other.p[1];
    this->p[2] += other.p[2];
    return *this;
  }

  triangle &operator-=(const triangle &other) {
    this->p[0] -= other.p[0];
    this->p[1] -= other.p[1];
    this->p[2] -= other.p[2];
    return *this;
  }

  triangle &operator*=(const triangle &other) {
    this->p[0] *= other.p[0];
    this->p[1] *= other.p[1];
    this->p[2] *= other.p[2];
    return *this;
  }

  triangle &operator/=(const triangle &other) {
    this->p[0] /= other.p[0];
    this->p[1] /= other.p[1];
    this->p[2] /= other.p[2];
    return *this;
  }

  std::vector<triangle> clip_against_plane(const vec3d &plane_p,
                                           vec3d plane_n) {
    std::vector<triangle> triangles;
    triangles.reserve(2);
    triangle temp;

    // Make sure plane normal is indeed normal
    plane_n.normalize();

    // Return signed shortest distance from point to plane, plane normal must be
    // normalised
    auto dist = [&](vec3d &p) {
      return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z -
              plane_n.dot_product(plane_p));
    };

    // Create two temporary storage arrays to classify points either side of
    // plane If distance sign is positive, point lies on "inside" of plane
    vec3d *inside_points[3];
    int nInsidePointCount = 0;
    vec3d *outside_points[3];
    int nOutsidePointCount = 0;

    vec2d *inside_tex[3];
    int nInsideTexCount = 0;
    vec2d *outside_tex[3];
    int nOutsideTexCount = 0;

    // Get signed distance of each point in triangle to plane
    float d0 = dist(this->p[0]);
    float d1 = dist(this->p[1]);
    float d2 = dist(this->p[2]);

    if (d0 >= 0) {
      inside_points[nInsidePointCount++] = &this->p[0];
      inside_tex[nInsideTexCount++] = &this->t[0];
    } else {
      outside_points[nOutsidePointCount++] = &this->p[0];
      outside_tex[nOutsideTexCount++] = &this->t[0];
    }
    if (d1 >= 0) {
      inside_points[nInsidePointCount++] = &this->p[1];
      inside_tex[nInsideTexCount++] = &this->t[1];
    } else {
      outside_points[nOutsidePointCount++] = &this->p[1];
      outside_tex[nOutsideTexCount++] = &this->t[1];
    }
    if (d2 >= 0) {
      inside_points[nInsidePointCount++] = &this->p[2];
      inside_tex[nInsideTexCount++] = &this->t[2];
    } else {
      outside_points[nOutsidePointCount++] = &this->p[2];
      outside_tex[nOutsideTexCount++] = &this->t[2];
    }

    // Now classify triangle points, and break the input triangle into
    // smaller output triangles if required. There are four possible
    // outcomes...

    if (nInsidePointCount == 3) {
      // All points lie on the inside of plane, so do nothing
      // and allow the triangle to simply pass through
      triangles.emplace_back(*this);
    }

    if (nInsidePointCount == 1 && nOutsidePointCount == 2) {
      // Triangle should be clipped. As two points lie outside
      // the plane, the triangle simply becomes a smaller triangle

      // Copy appearance info to new triangle
      temp.col = this->col;

      // The inside point is valid, so keep that...
      temp.p[0] = *inside_points[0];
      temp.t[0] = *inside_tex[0];

      // but the two new points are at the locations where the
      // original sides of the triangle (lines) intersect with the plane
      float t;
      temp.p[1] = vector_intersect_plane(plane_p, plane_n, *inside_points[0],
                                         *outside_points[0], t);
      temp.t[1].u =
          t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
      temp.t[1].v =
          t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;

      temp.p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[0],
                                         *outside_points[1], t);
      temp.t[2].u =
          t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
      temp.t[2].v =
          t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;

      triangles.emplace_back(temp);
    }

    if (nInsidePointCount == 2 && nOutsidePointCount == 1) {
      // Triangle should be clipped. As two points lie inside the plane,
      // the clipped triangle becomes a "quad". Fortunately, we can
      // represent a quad with two new triangles

      // Copy appearance info to new triangles
      temp.col = this->col;

      // The first triangle consists of the two inside points and a new
      // point determined by the location where one side of the triangle
      // intersects with the plane
      temp.p[0] = *inside_points[0];
      temp.p[1] = *inside_points[1];

      temp.t[0] = *inside_tex[0];
      temp.t[1] = *inside_tex[1];

      float t;

      temp.p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[0],
                                         *outside_points[0], t);
      temp.t[2].u =
          t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
      temp.t[2].v =
          t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;

      triangles.emplace_back(temp);

      // The second triangle is composed of one of he inside points, a
      // new point determined by the intersection of the other side of the
      // triangle and the plane, and the newly created point above
      temp.p[0] = *inside_points[1];
      temp.p[1] = temp.p[2];

      temp.t[0] = *inside_tex[1];
      temp.t[1] = temp.t[2];

      temp.p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[1],
                                         *outside_points[0], t);
      temp.t[2].u =
          t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
      temp.t[2].v =
          t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;

      triangles.emplace_back(temp);
    }
    return triangles;
  }
};

struct mesh {
  std::vector<triangle> tris;

  bool loadObj(std::string sPath, bool bHasTexture = false) {
    std::ifstream f(sPath);
    if (!f.is_open() || !f.good())
      return false;

    tris.reserve(RESERVE_TRIS);

    std::vector<vec3d> verts;
    verts.reserve(RESERVE_VERTS);

    std::vector<vec2d> texs;
    verts.reserve(RESERVE_VERTS);

    while (!f.eof()) {
      char line[128];
      f.getline(line, 128);

      std::strstream s;
      s << line;

      char junk;

      if (line[0] == 'v') {
        if (line[1] == 't') {
          vec2d v;
          s >> junk >> junk >> v.u >> v.v;
          texs.emplace_back(v);
        } else {
          vec3d v;
          s >> junk >> v.x >> v.y >> v.z;
          verts.emplace_back(v);
        }
      } else if (line[0] == 'f' && !bHasTexture) {
        int f[3];
        s >> junk >> f[0] >> f[1] >> f[2];
        tris.emplace_back(verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1]);
      } else if (line[0] == 'f' && bHasTexture) {
        s >> junk;

        std::string tokens[6];
        int nTokenCount = -1;

        while (!s.eof()) {
          char c = s.get();
          if (c == ' ' || c == '/')
            nTokenCount++;
          else
            tokens[nTokenCount].append(1, c);
        }

        tokens[nTokenCount].pop_back();

        tris.emplace_back(verts[stoi(tokens[0]) - 1],
                          verts[stoi(tokens[2]) - 1],
                          verts[stoi(tokens[4]) - 1], texs[stoi(tokens[1]) - 1],
                          texs[stoi(tokens[3]) - 1], texs[stoi(tokens[5]) - 1]);
      }
    }

    return true;
  }
};

olc::Pixel GetColour(float lum) {
  int nValue = (int)(std::max(lum, 0.20f) * 255.0f);
  return olc::Pixel(nValue, nValue, nValue);
}

class olcGameEngine : public olc::PixelGameEngine {
public:
  olcGameEngine() { sAppName = "3D demo"; }

private:
  mesh mMesh;
  mat4 matProj;

  vec3d vCamera;
  vec3d vLookDir;

  float fYaw;
  // float fTheta;
  float *pDepthBuffer;

  // olc::Sprite *sprTex;

public:
  bool OnUserCreate() override {
    // pDepthBuffer = new float[ScreenWidth() * ScreenHeight()];

    /*
        // SOUTH
        {{0.0f, 0.0f, 0.0f, 1.0f},
         {0.0f, 1.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{0.0f, 0.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 0.0f, 1.0f},
         {1.0f, 0.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},

        // EAST
        {{1.0f, 0.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{1.0f, 0.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 1.0f, 1.0f},
         {1.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},

        // NORTH
        {{1.0f, 0.0f, 1.0f, 1.0f},
         {1.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{1.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},

        // WEST
        {{0.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 1.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{0.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 1.0f, 0.0f, 1.0f},
         {0.0f, 0.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},

        // TOP
        {{0.0f, 1.0f, 0.0f, 1.0f},
         {0.0f, 1.0f, 1.0f, 1.0f},
         {1.0f, 1.0f, 1.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{0.0f, 1.0f, 0.0f, 1.0f},
         {1.0f, 1.0f, 1.0f, 1.0f},
         {1.0f, 1.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},

        // BOTTOM
        {{1.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 0.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {0.0f, 0.0f},
         {1.0f, 0.0f}},
        {{1.0f, 0.0f, 1.0f, 1.0f},
         {0.0f, 0.0f, 0.0f, 1.0f},
         {1.0f, 0.0f, 0.0f, 1.0f},
         {0.0f, 1.0f},
         {1.0f, 0.0f},
         {1.0f, 1.0f}},
    */

    mMesh.loadObj("mountains.obj");

    // sprTex = new olc::Sprite("jario.png");

    matProj.projection(-90.0f,
                       static_cast<float>(ScreenHeight()) /
                           static_cast<float>(ScreenWidth()),
                       NEAR_PLANE, 1000.0f);

    return true;
  }

  bool OnUserUpdate(float fElapsedTime) override {
    if (GetKey(olc::Key::UP).bHeld)
      vCamera.y += 8.0f * fElapsedTime;
    if (GetKey(olc::Key::DOWN).bHeld)
      vCamera.y -= 8.0f * fElapsedTime;
    if (GetKey(olc::Key::LEFT).bHeld)
      vCamera.x += 8.0f * fElapsedTime;
    if (GetKey(olc::Key::RIGHT).bHeld)
      vCamera.x -= 8.0f * fElapsedTime;

    vec3d vForward = vLookDir * (8.0f * fElapsedTime);
    if (GetKey(olc::Key::W).bHeld)
      vCamera += vForward;
    if (GetKey(olc::Key::S).bHeld)
      vCamera -= vForward;
    if (GetKey(olc::Key::A).bHeld)
      fYaw -= 2.0f * fElapsedTime;
    if (GetKey(olc::Key::D).bHeld)
      fYaw += 2.0f * fElapsedTime;

    // mat4 matRotZ, matRotX;
    // fTheta += fElapsedTime;
    // matRotZ.rotation_z(fTheta / 2.0f);
    // matRotX.rotation_x(fTheta);

    // Make translation matrix
    mat4 matTrans;
    matTrans.translation(0.0f, 0.0f, 8.0f);

    // Make world matrix
    mat4 matWorld;
    matWorld.identity();
    // matWorld = matRotZ * matRotX;
    matWorld = matWorld * matTrans;

    vec3d vUp = {0.0f, 1.0f, 0.0f};
    vec3d vTarget = {0.0f, 0.0f, 1.0f};
    mat4 matCameraRot;
    matCameraRot.rotation_y(fYaw);
    vLookDir = vTarget * matCameraRot;
    vTarget = vCamera + vLookDir;

    mat4 matCamera;
    matCamera.point_at(vCamera, vTarget, vUp);
    matCamera.quick_inverse();

    // Make vector of trinagles to raster
    std::vector<triangle> vecTriangleToRaster;
    vecTriangleToRaster.reserve(RESERVE_RAST_TRIS);

    // Draw triangles
    for (triangle &tri : mMesh.tris) {
      triangle triProjected, triTransformed, triViewed;

      // Rotate triangles
      triTransformed.p[0] = tri.p[0] * matWorld;
      triTransformed.p[1] = tri.p[1] * matWorld;
      triTransformed.p[2] = tri.p[2] * matWorld;
      triTransformed.t[0] = tri.t[0];
      triTransformed.t[1] = tri.t[1];
      triTransformed.t[2] = tri.t[2];

      // Calculate triangle normal
      vec3d normal, line1, line2;

      // Get either side of triangle
      line1 = triTransformed.p[1] - triTransformed.p[0];
      line2 = triTransformed.p[2] - triTransformed.p[0];

      // Take cross product to get normal
      normal = line1.cross_product(line2);
      normal.normalize();

      // Get ray from triangle to camera
      vec3d vCameraRay = triTransformed.p[0] - vCamera;

      // If ray is alligned with normal then triangle is visible
      if (normal.dot_product(vCameraRay) < 0.0f) {
        // Move triangles from world space to view space
        triViewed.p[0] = triTransformed.p[0] * matCamera;
        triViewed.p[1] = triTransformed.p[1] * matCamera;
        triViewed.p[2] = triTransformed.p[2] * matCamera;
        triViewed.t[0] = triTransformed.t[0];
        triViewed.t[1] = triTransformed.t[1];
        triViewed.t[2] = triTransformed.t[2];

        // Clip triViewed against near plane
        std::vector<triangle> clipped = triViewed.clip_against_plane(
            {0.0f, 0.0f, NEAR_PLANE}, {0.0f, 0.0f, 1.0f});

        olc::Pixel c;
        if (clipped.size()) {
          // Set and normalize light direction
          vec3d light_direction = {0.0f, 1.0f, -1.0f};
          light_direction.normalize();

          // How alligned are trangle surface normal and light direction
          float dp = std::max(0.1f, light_direction.dot_product(normal));

          c = GetColour(dp);
        }

        for (size_t i = 0; i < clipped.size(); i++) {
          // Project trinagles from 3D to 2D
          triProjected.p[0] = clipped[i].p[0] * matProj;
          triProjected.p[1] = clipped[i].p[1] * matProj;
          triProjected.p[2] = clipped[i].p[2] * matProj;
          triProjected.col = c;
          triProjected.t[0] = clipped[i].t[0];
          triProjected.t[1] = clipped[i].t[1];
          triProjected.t[2] = clipped[i].t[2];

          triProjected.t[0] /= triProjected.p[0].w;
          triProjected.t[1] /= triProjected.p[1].w;
          triProjected.t[2] /= triProjected.p[2].w;

          triProjected.t[0].w = 1.0f / triProjected.p[0].w;
          triProjected.t[1].w = 1.0f / triProjected.p[1].w;
          triProjected.t[2].w = 1.0f / triProjected.p[2].w;

          triProjected.p[0] /= triProjected.p[0].w;
          triProjected.p[1] /= triProjected.p[1].w;
          triProjected.p[2] /= triProjected.p[2].w;

          // Scale into view
          vec3d vOffsetView = {1.0f, 1.0f, 0.0f};
          triProjected.p[0] += vOffsetView;
          triProjected.p[1] += vOffsetView;
          triProjected.p[2] += vOffsetView;

          triProjected.p[0].x *= 0.5f * static_cast<float>(ScreenWidth());
          triProjected.p[1].x *= 0.5f * static_cast<float>(ScreenWidth());
          triProjected.p[2].x *= 0.5f * static_cast<float>(ScreenWidth());
          triProjected.p[0].y *= 0.5f * static_cast<float>(ScreenHeight());
          triProjected.p[1].y *= 0.5f * static_cast<float>(ScreenHeight());
          triProjected.p[2].y *= 0.5f * static_cast<float>(ScreenHeight());

          vecTriangleToRaster.emplace_back(triProjected);
        }
      }
    }

    // Sort triangles by depth (from back to front)
    std::sort(vecTriangleToRaster.begin(), vecTriangleToRaster.end(),
              [](const triangle &t1, const triangle &t2) {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
                return z1 > z2;
              });

    FillRect(0, 0, ScreenWidth(), ScreenHeight(), olc::BLACK);
    /*
    for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++)
      pDepthBuffer[i] = 0.0f;
    */

    for (triangle &triToRaster : vecTriangleToRaster) {
      // Clip triangles against all four screen edges, this could yield
      // a bunch of triangles, so create a queue that we traverse to
      //  ensure we only test new triangles generated against planes
      std::vector<triangle> clipped;
      clipped.reserve(2);
      std::list<triangle> listTriangles;

      // Add initial triangle
      listTriangles.emplace_back(triToRaster);
      int nNewTriangles = 1;

      for (int p = 0; p < 4; p++) {
        while (nNewTriangles > 0) {
          // Take triangle from front of queue
          triangle test = listTriangles.front();
          listTriangles.pop_front();
          nNewTriangles--;

          // Clip it against a plane. We only need to test each
          // subsequent plane, against subsequent new triangles
          // as all triangles after a plane clip are guaranteed
          // to lie on the inside of the plane. I like how this
          // comment is almost completely and utterly justified
          switch (p) {
          case 0:
            clipped =
                test.clip_against_plane({0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
            break;
          case 1:
            clipped = test.clip_against_plane(
                {0.0f, static_cast<float>(ScreenHeight()) - 1, 0.0f},
                {0.0f, -1.0f, 0.0f});
            break;
          case 2:
            clipped =
                test.clip_against_plane({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
            break;
          case 3:
            clipped = test.clip_against_plane(
                {static_cast<float>(ScreenWidth()) - 1, 0.0f, 0.0f},
                {-1.0f, 0.0f, 0.0f});
            break;
          }

          // Clipping may yield a variable number of triangles, so
          // add these new ones to the back of the queue for subsequent
          // clipping against next planes
          for (size_t w = 0; w < clipped.size(); w++)
            listTriangles.emplace_back(clipped[w]);
        }
        nNewTriangles = listTriangles.size();
      }

      // Draw the transformed, viewed, clipped, projected, sorted, clipped
      // triangles
      for (auto &t : listTriangles) {
        /*
        textured_triangle(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
                          t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
                          t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w,
                          sprTex);
        */
        FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y,
                     t.col);
#ifdef DEBUG
        DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y,
                     olc::WHITE);
#endif
      }
    }

    return true;
  }

  void textured_triangle(int x1, int y1, float u1, float v1, float w1, int x2,
                         int y2, float u2, float v2, float w2, int x3, int y3,
                         float u3, float v3, float w3,
                         const olc::Sprite *const tex) {
    if (y2 < y1) {
      std::swap(y1, y2);
      std::swap(x1, x2);
      std::swap(u1, u2);
      std::swap(v1, v2);
      std::swap(w1, w2);
    }

    if (y3 < y1) {
      std::swap(y1, y3);
      std::swap(x1, x3);
      std::swap(u1, u3);
      std::swap(v1, v3);
      std::swap(w1, w3);
    }

    if (y3 < y2) {
      std::swap(y2, y3);
      std::swap(x2, x3);
      std::swap(u2, u3);
      std::swap(v2, v3);
      std::swap(w2, w3);
    }

    int dy1 = y2 - y1;
    int dx1 = x2 - x1;
    float du1 = u2 - u1;
    float dv1 = v2 - v1;
    float dw1 = w2 - w1;

    int dy2 = y3 - y1;
    int dx2 = x3 - x1;
    float du2 = u3 - u1;
    float dv2 = v3 - v1;
    float dw2 = w3 - w1;

    float tex_u, tex_v, tex_w;

    float dax_step = 0.0f;
    float dbx_step = 0.0f;
    float du1_step = 0.0f;
    float dv1_step = 0.0f;
    float dw1_step = 0.0f;
    float du2_step = 0.0f;
    float dv2_step = 0.0f;
    float dw2_step = 0.0f;

    if (dy1) {
      dax_step = dx1 / static_cast<float>(std::abs(dy1));
      du1_step = du1 / static_cast<float>(std::abs(dy1));
      dv1_step = dv1 / static_cast<float>(std::abs(dy1));
      dw1_step = dw1 / static_cast<float>(std::abs(dy1));
    }
    if (dy2) {
      dbx_step = dx2 / static_cast<float>(std::abs(dy2));
      du2_step = du2 / static_cast<float>(std::abs(dy2));
      dv2_step = dv2 / static_cast<float>(std::abs(dy2));
      dw2_step = dw2 / static_cast<float>(std::abs(dy2));
    }

    if (dy1) {
      for (int i = y1; i <= y2; i++) {
        int ax = x1 + static_cast<float>(i - y1) * dax_step;
        int bx = x1 + static_cast<float>(i - y1) * dbx_step;

        float tex_su = u1 + static_cast<float>(i - y1) * du1_step;
        float tex_sv = v1 + static_cast<float>(i - y1) * dv1_step;
        float tex_sw = w1 + static_cast<float>(i - y1) * dw1_step;

        float tex_eu = u1 + static_cast<float>(i - y1) * du2_step;
        float tex_ev = v1 + static_cast<float>(i - y1) * dv2_step;
        float tex_ew = w1 + static_cast<float>(i - y1) * dw2_step;

        if (ax > bx) {
          std::swap(ax, bx);
          std::swap(tex_su, tex_eu);
          std::swap(tex_sv, tex_ev);
          std::swap(tex_sw, tex_ew);
        }

        float tstep = 1.0f / static_cast<float>(bx - ax);
        float t = 0.0f;

        for (int j = ax; j < bx; j++) {
          tex_u = (1.0f - t) * tex_su + t * tex_eu;
          tex_v = (1.0f - t) * tex_sv + t * tex_ev;
          tex_w = (1.0f - t) * tex_sw + t * tex_ew;

          if (tex_w > pDepthBuffer[i * ScreenWidth() + j]) {
            Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
            pDepthBuffer[i * ScreenWidth() + j] = tex_w;
          }

          t += tstep;
        }
      }
    }

    dy1 = y3 - y2;
    dx1 = x3 - x2;
    du1 = u3 - u2;
    dv1 = v3 - v2;
    dw1 = w3 - w2;

    du1_step = 0.0f;
    dv1_step = 0.0f;

    if (dy1) {
      dax_step = dx1 / static_cast<float>(std::abs(dy1));
      du1_step = du1 / static_cast<float>(std::abs(dy1));
      dv1_step = dv1 / static_cast<float>(std::abs(dy1));
      dw1_step = dw1 / static_cast<float>(std::abs(dy1));
    }
    if (dy2) {
      dbx_step = dx2 / static_cast<float>(std::abs(dy2));
    }

    if (dy1) {
      for (int i = y2; i <= y3; i++) {
        int ax = x2 + static_cast<float>(i - y2) * dax_step;
        int bx = x1 + static_cast<float>(i - y1) * dbx_step;

        float tex_su = u2 + static_cast<float>(i - y2) * du1_step;
        float tex_sv = v2 + static_cast<float>(i - y2) * dv1_step;
        float tex_sw = w2 + static_cast<float>(i - y2) * dw1_step;

        float tex_eu = u1 + static_cast<float>(i - y1) * du2_step;
        float tex_ev = v1 + static_cast<float>(i - y1) * dv2_step;
        float tex_ew = w1 + static_cast<float>(i - y1) * dw2_step;

        if (ax > bx) {
          std::swap(ax, bx);
          std::swap(tex_su, tex_eu);
          std::swap(tex_sv, tex_ev);
          std::swap(tex_sw, tex_ew);
        }

        float tstep = 1.0f / static_cast<float>(bx - ax);
        float t = 0.0f;

        for (int j = ax; j < bx; j++) {
          tex_u = (1.0f - t) * tex_su + t * tex_eu;
          tex_v = (1.0f - t) * tex_sv + t * tex_ev;
          tex_w = (1.0f - t) * tex_sw + t * tex_ew;

          if (tex_w > pDepthBuffer[i * ScreenWidth() + j]) {
            Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
            pDepthBuffer[i * ScreenWidth() + j] = tex_w;
          }

          t += tstep;
        }
      }
    }
  }
};

int main() {
  olcGameEngine demo;
  if (demo.Construct(256, 240, 4, 4))
    demo.Start();
  else
    std::cerr << "Could not construct" << std::endl;

  return EXIT_SUCCESS;
}
