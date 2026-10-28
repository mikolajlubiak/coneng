# ConEng
3D (no longer CONsole) ENGine

A software 3D renderer in C++: the whole triangle pipeline runs on the CPU, and the finished frame is handed to [olcPixelGameEngine](https://github.com/OneLoneCoder/olcPixelGameEngine) only for display. It started as a console renderer, hence the name.

Demo: https://video.infosec.exchange/w/p/f3rbPWUWt5zTMZsLAJ3SMP

## Pipeline

- Model → world → view → projection transforms on operator-overloaded `vec3d` / `mat4` types, with a look-at camera (fast inverse) and a perspective projection
- Back-face culling and flat Lambert shading from a directional light
- Clipping against the near plane in view space, then against the four screen edges
- Perspective-correct texture mapping (u/w, v/w and 1/w interpolation) with a 1/w depth buffer for textured triangles
- Painter's-algorithm sorting for flat-shaded meshes
- OBJ loading (vertex positions, texture coordinates and triangle faces)
- Several meshes with independent world transforms in one scene: a spinning textured cube flying over a flat-shaded terrain
- Free camera with yaw, pitch and mouse look

## Controls

- `W` / `S`: move forward / backward
- `A` / `D`: turn left / right
- `Page Up` / `Page Down`: look up / down
- Arrow keys: move the camera along the world X and Y axes
- Mouse: look around
- `Esc`: quit

## Build

Linux:

- Install a C++ compiler, make, and the X11, OpenGL and libpng headers:
  - Ubuntu: `sudo apt install g++ make libx11-dev libgl-dev libpng-dev`
  - Fedora: `sudo dnf install gcc-c++ make libX11-devel mesa-libGL-devel libpng-devel`
- `make release` (or `make debug`)
- `./release` (run it from the repository root so the models and texture are found)
