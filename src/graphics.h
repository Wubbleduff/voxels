
#pragma once

#include "common.h"
#include "math.h"

f32 get_screen_aspect_ratio();

void draw_mesh(
        u32 num_verts,
        f32* vx,
        f32* vy,
        f32* vz,
        f32* nx,
        f32* ny,
        f32* nz,
        u32 num_indices,
        u16* indices,
        const mtx4x4* mvp);

void debug_draw_line(const mtx4x4* camera_and_clip_mtx, v3 a, v3 b, v3 c);
void debug_draw_sphere(const mtx4x4* proj_mtx, const mtx4x4* cam_mtx, v3 pos, f32 r, v3 c);


