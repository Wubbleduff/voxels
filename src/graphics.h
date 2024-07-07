
#pragma once

#include "common.h"
#include "math.h"

#define MESH_64_MAX_VERTICES (64)
#define MESH_64_MAX_INDICES (2*64)
struct Mesh_64
{
    u32_m num_vertices;
    f32_m vx[MESH_64_MAX_VERTICES];
    f32_m vy[MESH_64_MAX_VERTICES];
    f32_m vz[MESH_64_MAX_VERTICES];
    f32_m nx[MESH_64_MAX_VERTICES];
    f32_m ny[MESH_64_MAX_VERTICES];
    f32_m nz[MESH_64_MAX_VERTICES];

    u32_m num_indices;
    u32_m indices[MESH_64_MAX_INDICES];
};

#define MESH_1M_MAX_VERTICES (1024*1024)
#define MESH_1M_MAX_INDICES (2*1024*1024)
struct Mesh_1M
{
    u32_m num_vertices;
    f32_m vx[MESH_1M_MAX_VERTICES];
    f32_m vy[MESH_1M_MAX_VERTICES];
    f32_m vz[MESH_1M_MAX_VERTICES];
    f32_m nx[MESH_1M_MAX_VERTICES];
    f32_m ny[MESH_1M_MAX_VERTICES];
    f32_m nz[MESH_1M_MAX_VERTICES];

    u32_m num_indices;
    u32_m indices[MESH_1M_MAX_INDICES];
};

struct FrameBuffer
{
    u32_m width;
    u32_m height;
    // Pixels stored as 0xRRGGBBAA
    u32_m* buf;
};

f32 get_screen_aspect_ratio();

void draw_Mesh_64(const struct Mesh_64* mesh, const mtx4x4* camera_and_clip_mtx, u64 num, f32* pos_x, f32* pos_y, f32* pos_z);
void draw_Mesh_1M(const struct Mesh_1M* mesh, const mtx4x4* camera_and_clip_mtx, u64 num, f32* pos_x, f32* pos_y, f32* pos_z);

void debug_draw_line(const mtx4x4* camera_and_clip_mtx, v3 a, v3 b, v3 c);


