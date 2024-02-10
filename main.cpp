#include "olcConsoleGameEngineSDL.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <strstream>
#include <vector>

#define RESERVE_VERTS 64
#define RESERVE_TRIS RESERVE_VERTS * 2
#define RESERVE_RAST_TRIS RESERVE_TRIS / 4

struct vec3d {
  float x, y, z = 0.0f;
  float w = 1.0f;

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

  float dot_product(const vec3d &other) {
    return this->x * other.x + this->y * other.y + this->z * other.z;
  }

  float length() { return sqrtf(this->dot_product(*this)); }

  void normalize() { *this /= this->length(); }

  vec3d cross_product(const vec3d &other) {
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

vec3d operator*(vec3d lhs, float rhs) {
  lhs.x *= rhs;
  lhs.y *= rhs;
  lhs.z *= rhs;
  return lhs;
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

    mat4 matrix;
    matrix.m[0][0] = newRight.x;
    matrix.m[0][1] = newRight.y;
    matrix.m[0][2] = newRight.z;
    matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = up.x;
    matrix.m[1][1] = up.y;
    matrix.m[1][2] = up.z;
    matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = target.x;
    matrix.m[2][1] = target.y;
    matrix.m[2][2] = target.z;
    matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;
    matrix.m[3][1] = pos.y;
    matrix.m[3][2] = pos.z;
    matrix.m[3][3] = 1.0f;
    *this = matrix;
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

  wchar_t sym = PIXEL_SOLID;
  short col = FG_WHITE;

#ifdef DEBUG
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
};

struct mesh {
  std::vector<triangle> tris;

  bool loadObj(std::string sPath) {
    std::ifstream f(sPath);
    if (!f.is_open() || !f.good())
      return false;

    tris.reserve(RESERVE_TRIS);

    std::vector<vec3d> verts;
    verts.reserve(RESERVE_VERTS);

    while (!f.eof()) {
      char line[128];
      f.getline(line, 128);

      std::strstream s;
      s << line;

      char junk;

      if (line[0] == 'v') {
        vec3d v;
        s >> junk >> v.x >> v.y >> v.z;
        verts.emplace_back(v);
      } else if (line[0] == 'f') {
        int f[3];
        s >> junk >> f[0] >> f[1] >> f[2];
        tris.push_back({{verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1]}});
      }
    }

    return true;
  }
};

CHAR_INFO GetColour(float lum) {
  short bg_col, fg_col;
  wchar_t sym;
  int pixel_bw = (int)(13.0f * lum);
  switch (pixel_bw) {
  case 0:
    bg_col = BG_BLACK;
    fg_col = FG_BLACK;
    sym = PIXEL_SOLID;
    break;

  case 1:
    bg_col = BG_BLACK;
    fg_col = FG_DARK_GREY;
    sym = PIXEL_QUARTER;
    break;
  case 2:
    bg_col = BG_BLACK;
    fg_col = FG_DARK_GREY;
    sym = PIXEL_HALF;
    break;
  case 3:
    bg_col = BG_BLACK;
    fg_col = FG_DARK_GREY;
    sym = PIXEL_THREEQUARTERS;
    break;
  case 4:
    bg_col = BG_BLACK;
    fg_col = FG_DARK_GREY;
    sym = PIXEL_SOLID;
    break;

  case 5:
    bg_col = BG_DARK_GREY;
    fg_col = FG_GREY;
    sym = PIXEL_QUARTER;
    break;
  case 6:
    bg_col = BG_DARK_GREY;
    fg_col = FG_GREY;
    sym = PIXEL_HALF;
    break;
  case 7:
    bg_col = BG_DARK_GREY;
    fg_col = FG_GREY;
    sym = PIXEL_THREEQUARTERS;
    break;
  case 8:
    bg_col = BG_DARK_GREY;
    fg_col = FG_GREY;
    sym = PIXEL_SOLID;
    break;

  case 9:
    bg_col = BG_GREY;
    fg_col = FG_WHITE;
    sym = PIXEL_QUARTER;
    break;
  case 10:
    bg_col = BG_GREY;
    fg_col = FG_WHITE;
    sym = PIXEL_HALF;
    break;
  case 11:
    bg_col = BG_GREY;
    fg_col = FG_WHITE;
    sym = PIXEL_THREEQUARTERS;
    break;
  case 12:
    bg_col = BG_GREY;
    fg_col = FG_WHITE;
    sym = PIXEL_SOLID;
    break;
  default:
    bg_col = BG_BLACK;
    fg_col = FG_BLACK;
    sym = PIXEL_SOLID;
  }

  CHAR_INFO c;
  c.colour = bg_col | fg_col;
  c.glyph = sym;
  return c;
}

class olcGameEngine : public olcConsoleGameEngine {
public:
  olcGameEngine() { m_sAppName = L"3D demo"; }

private:
  mesh mMesh;
  mat4 matProj;

  vec3d vCamera;
  vec3d vLookDir;

  float fYaw;

public:
  bool OnUserCreate() override {
    mMesh.loadObj("axis.obj");

    matProj.projection(-90.0f,
                       static_cast<float>(ScreenHeight()) /
                           static_cast<float>(ScreenWidth()),
                       0.1f, 1000.0f);

    return true;
  }

  bool OnUserUpdate(float fElapsedTime) override {
    if (GetKey(VK_UP).bHeld)
      vCamera.y += 8.0f * fElapsedTime;
    if (GetKey(VK_DOWN).bHeld)
      vCamera.y -= 8.0f * fElapsedTime;
    if (GetKey(VK_LEFT).bHeld)
      vCamera.x += 8.0f * fElapsedTime;
    if (GetKey(VK_RIGHT).bHeld)
      vCamera.x -= 8.0f * fElapsedTime;

    vec3d vForward = vLookDir * (8.0f * fElapsedTime);
    if (GetKey(L'W').bHeld)
      vCamera += vForward;
    if (GetKey(L'S').bHeld)
      vCamera -= vForward;
    if (GetKey(L'A').bHeld)
      fYaw -= 2.0f * fElapsedTime;
    if (GetKey(L'D').bHeld)
      fYaw += 2.0f * fElapsedTime;

    Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

    // Make translation matrix
    mat4 matTrans;
    matTrans.translation(0.0f, 0.0f, 16.0f);

    // Make world matrix
    mat4 matWorld;
    matWorld.identity();
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
        // Set and normalize light direction
        vec3d light_direction = {0.0f, 1.0f, -1.0f};
        light_direction.normalize();

        // How alligned are trangle surface normal and light direction
        float dp = std::max(0.1f, light_direction.dot_product(normal));

        CHAR_INFO c = GetColour(dp);

        // Move triangles from world space to view space
        triViewed.p[0] = triTransformed.p[0] * matCamera;
        triViewed.p[1] = triTransformed.p[1] * matCamera;
        triViewed.p[2] = triTransformed.p[2] * matCamera;

        // Project trinagles from 3D to 2D
        triProjected.p[0] = triViewed.p[0] * matProj;
        triProjected.p[1] = triViewed.p[1] * matProj;
        triProjected.p[2] = triViewed.p[2] * matProj;
        triProjected.col = c.colour;
        triProjected.sym = c.glyph;

        triProjected.p[0] /= triProjected.p[0].w;
        triProjected.p[1] /= triProjected.p[1].w;
        triProjected.p[2] /= triProjected.p[2].w;

        // Scale into view
        vec3d vOffsetView = {1.0f, 1.0f, 0.0f};
        triProjected.p[0] += vOffsetView;
        triProjected.p[1] += vOffsetView;
        triProjected.p[2] += vOffsetView;

        triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
        triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
        triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
        triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
        triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
        triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

        vecTriangleToRaster.emplace_back(triProjected);
      }
    }

    // Sort triangles by depth (from back to front)
    std::sort(vecTriangleToRaster.begin(), vecTriangleToRaster.end(),
              [](const triangle &t1, const triangle &t2) {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
                return z1 > z2;
              });

    for (triangle &triProjected : vecTriangleToRaster) {
      FillTriangle(triProjected.p[0].x, triProjected.p[0].y,
                   triProjected.p[1].x, triProjected.p[1].y,
                   triProjected.p[2].x, triProjected.p[2].y, triProjected.sym,
                   triProjected.col);
#ifdef DEBUG
      DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
                   triProjected.p[1].x, triProjected.p[1].y,
                   triProjected.p[2].x, triProjected.p[2].y, PIXEL_SOLID,
                   FG_WHITE);
#endif
    }

    return true;
  }
};

int main() {
  olcGameEngine demo;
  if (demo.ConstructConsole(256, 240, 4, 4))
    demo.Start();
  else
    std::cerr << "Could not construct console" << std::endl;

  return 0;
}
