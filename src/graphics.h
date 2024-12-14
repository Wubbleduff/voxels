
#pragma once

#include "common.h"
#include "math.h"

f32 get_screen_aspect_ratio();

void draw_mesh(
        const u32 num_verts,
        const f32* vx,
        const f32* vy,
        const f32* vz,
        const f32* nx,
        const f32* ny,
        const f32* nz,
        const u32 num_indices,
        const u16* indices,
        const mtx4x4* mvp);

void debug_draw_line(const mtx4x4* camera_and_clip_mtx, v3 a, v3 b, v3 c);
void debug_draw_sphere(const mtx4x4* proj_mtx, const mtx4x4* cam_mtx, v3 pos, const f32 r, v3 c);


