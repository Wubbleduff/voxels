
#pragma once

#include "common.h"
#include "math.h"

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

void draw_imm_Mesh_1M(const struct Mesh_1M* p, const mtx4x4* transform);

void debug_draw_imm_line(const mtx4x4* transform, v3 a, v3 b, v3 c);


