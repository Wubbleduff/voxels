
#pragma once

#include "common.h"

struct MemoryArena;

#define TERRAIN_MAX_VERTICES (400000)
#define TERRAIN_MAX_INDICES (800000)
struct Terrain
{
    u32_m num_vertices;
    f32_m vx[TERRAIN_MAX_VERTICES];
    f32_m vy[TERRAIN_MAX_VERTICES];
    f32_m vz[TERRAIN_MAX_VERTICES];
    f32_m nx[TERRAIN_MAX_VERTICES];
    f32_m ny[TERRAIN_MAX_VERTICES];
    f32_m nz[TERRAIN_MAX_VERTICES];

    u32_m num_indices;
    u32_m indices[TERRAIN_MAX_INDICES];
};

INTERNAL void reset_terrain(struct Terrain* terrain)
{
    terrain->num_vertices = 0;
    terrain->num_indices = 0;
}

INTERNAL void copy_terrain(struct Terrain* dst, const struct Terrain* src)
{
    dst->num_vertices = src->num_vertices;
    ARRAY_COPY(dst->vx, src->vx, src->num_vertices);
    ARRAY_COPY(dst->vy, src->vy, src->num_vertices);
    ARRAY_COPY(dst->vz, src->vz, src->num_vertices);
    ARRAY_COPY(dst->nx, src->nx, src->num_vertices);
    ARRAY_COPY(dst->ny, src->ny, src->num_vertices);
    ARRAY_COPY(dst->nz, src->nz, src->num_vertices);
    dst->num_indices = src->num_indices;
    ARRAY_COPY(dst->indices, src->indices, src->num_indices);
}

void generate_terrain(
        u32_m* out_num_vertices,
        f32_m* out_vx,
        f32_m* out_vy,
        f32_m* out_vz,
        f32_m* out_nx,
        f32_m* out_ny,
        f32_m* out_nz,
        u32_m* out_num_indices,
        u32_m* out_indices,
        f32 cam_pos_x,
        f32 cam_pos_y,
        f32 cam_pos_z,
        u8 lod);
