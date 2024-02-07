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
  float x, y, z;
};

struct triangle {
  vec3d p[3];

  wchar_t sym;
  short col;

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

struct mat4 {
  float m[4][4] = {0.0f};
};

class olcGameEngine : public olcConsoleGameEngine {

public:
  olcGameEngine() { m_sAppName = L"3D demo"; }

private:
  mesh mMesh;
  mat4 matProj;

  vec3d vCamera;

  float fTheta;

  void MultiplyMatrixVector(vec3d &i, vec3d &o, mat4 &m) {
    o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
    o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
    o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
    float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

    if (w != 0.0f) {
      o.x = o.x / w;
      o.y = o.y / w;
      o.z = o.z / w;
    }
  }

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

public:
  bool OnUserCreate() override {
    /*
    meshCube.tris = {

        // SOUTH
        {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f},

        // EAST
        {1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f},
        {1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f},

        // NORTH
        {1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f},
        {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f},

        // WEST
        {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},

        // TOP
        {0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f},

        // BOTTOM
        {1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f},
        {1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f},

    };
    */

    mMesh.loadObj("ship.obj");

    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 90.0f;
    float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
    float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.141592f);

    matProj.m[0][0] = fAspectRatio * fFovRad;
    matProj.m[1][1] = fFovRad;
    matProj.m[2][2] = fFar / (fFar - fNear);
    matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
    matProj.m[2][3] = 1.0f;
    matProj.m[3][3] = 0.0f;

    return true;
  }

  bool OnUserUpdate(float fElapsedTime) override {
    Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

    mat4 matRotZ, matRotX;

    fTheta += 1.0f * fElapsedTime;

    // Rotation Z
    matRotZ.m[0][0] = cosf(fTheta);
    matRotZ.m[0][1] = sinf(fTheta);
    matRotZ.m[1][0] = -sinf(fTheta);
    matRotZ.m[1][1] = cosf(fTheta);
    matRotZ.m[2][2] = 1;
    matRotZ.m[3][3] = 1;

    // Rotation X
    matRotX.m[0][0] = 1;
    matRotX.m[1][1] = cosf(fTheta * 0.5f);
    matRotX.m[1][2] = sinf(fTheta * 0.5f);
    matRotX.m[2][1] = -sinf(fTheta * 0.5f);
    matRotX.m[2][2] = cosf(fTheta * 0.5f);
    matRotX.m[3][3] = 1;

    std::vector<triangle> vecTriangleToRaster;
    vecTriangleToRaster.reserve(RESERVE_RAST_TRIS);

    // Draw triangles
    for (triangle tri : mMesh.tris) {
      triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

      // Rotate triangles

      // Rotate Z
      MultiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
      MultiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
      MultiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

      // Rotate X
      MultiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
      MultiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
      MultiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);

      // Translate triangles
      triTranslated = triRotatedZX;
      triTranslated.p[0].z = triTranslated.p[0].z + 8.0f;
      triTranslated.p[1].z = triTranslated.p[1].z + 8.0f;
      triTranslated.p[2].z = triTranslated.p[2].z + 8.0f;

      // Check if side is visible
      vec3d normal, line1, line2;

      line1.x = triTranslated.p[1].x - triTranslated.p[0].x;
      line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
      line1.z = triTranslated.p[1].z - triTranslated.p[0].z;

      line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
      line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
      line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

      normal.x = line1.y * line2.z - line1.z * line2.y;
      normal.y = line1.z * line2.x - line1.x * line2.z;
      normal.z = line1.x * line2.y - line1.y * line2.x;

      float l = sqrt(pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2));
      normal.x /= l;
      normal.y /= l;
      normal.z /= l;

      if ((normal.x * (triTranslated.p[0].x - vCamera.x) +
           normal.y * (triTranslated.p[0].y - vCamera.y) +
           normal.z * (triTranslated.p[0].z - vCamera.z)) < 0.0f) {

        // Add basic lighting
        vec3d light_direction = {0.0f, 0.0f, -1.0f};

        float l = sqrt(pow(light_direction.x, 2) + pow(light_direction.y, 2) +
                       pow(light_direction.z, 2));
        light_direction.x /= l;
        light_direction.y /= l;
        light_direction.z /= l;

        float dp = normal.x * light_direction.x + normal.y * light_direction.y +
                   normal.z * light_direction.z;

        CHAR_INFO c = GetColour(dp);
        triTranslated.col = c.colour;
        triTranslated.sym = c.glyph;

        // Project trinagles from 3D to 2D
        MultiplyMatrixVector(triTranslated.p[0], triProjected.p[0], matProj);
        MultiplyMatrixVector(triTranslated.p[1], triProjected.p[1], matProj);
        MultiplyMatrixVector(triTranslated.p[2], triProjected.p[2], matProj);
        triProjected.col = triTranslated.col;
        triProjected.sym = triTranslated.sym;

        // Scale into view
        triProjected.p[0].x += 1.0f;
        triProjected.p[0].y += 1.0f;
        triProjected.p[1].x += 1.0f;
        triProjected.p[1].y += 1.0f;
        triProjected.p[2].x += 1.0f;
        triProjected.p[2].y += 1.0f;

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
