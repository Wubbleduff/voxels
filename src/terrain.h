
#pragma once

#include "common.h"

struct MemoryArena;

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
        f32 cam_pos_z);
