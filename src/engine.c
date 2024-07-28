
#include "common.h"
#include "platform.h"
#include "graphics.h"
#include "math.h"

#include "marching_cubes.h"


INTERNAL inline __m256 pnoise8_calc_gradient(__m256 x, __m256 y, __m256 z, __m256 vx, __m256 vy, __m256 vz)
{
    u32 num_sphere_points = 1 << 16;
    // Create noise
    __m256i i_noise = rand8_u32(_mm256_castps_si256(_mm256_xor_ps(_mm256_xor_ps(vx, vy), vz)));
    i_noise = _mm256_and_si256(i_noise, _mm256_set1_epi32(num_sphere_points - 1));
    // Use random number to index points on a sphere.
    const __m256 noise = _mm256_cvtepi32_ps(i_noise);
    const __m256 u = _mm256_fmsub_ps(_mm256_set1_ps(2.0f / (f32)(num_sphere_points - 1)), noise, _mm256_set1_ps(1.0f));
    const __m256 t = _mm256_mul_ps(_mm256_set1_ps(10.166640738f), noise);
    const __m256 up = _mm256_sqrt_ps(
            _mm256_max_ps(
                _mm256_set1_ps(0.0f),
                _mm256_fnmadd_ps(u, u, _mm256_set1_ps(1.0f))
                )
            );
    const __m256 gx = _mm256_mul_ps(up, approx_cos8(t));
    const __m256 gy = _mm256_mul_ps(up, approx_sin8(t));
    const __m256 gz = u;

    __m256 dx = _mm256_sub_ps(x, vx);
    __m256 dy = _mm256_sub_ps(y, vy);
    __m256 dz = _mm256_sub_ps(z, vz);
    const __m256 result = _mm256_fmadd_ps(gx, dx, _mm256_fmadd_ps(gy, dy, _mm256_mul_ps(gz, dz)));
    return result;
}

INTERNAL inline __m256 pnoise8(const __m256 x, const __m256 y, const __m256 z)
{
    __m256 x0 = _mm256_round_ps(x, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
    __m256 x1 = _mm256_round_ps(_mm256_add_ps(x, _mm256_set1_ps(1.0f)), (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
    __m256 y0 = _mm256_round_ps(y, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
    __m256 y1 = _mm256_round_ps(_mm256_add_ps(y, _mm256_set1_ps(1.0f)), (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
    __m256 z0 = _mm256_round_ps(z, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
    __m256 z1 = _mm256_round_ps(_mm256_add_ps(z, _mm256_set1_ps(1.0f)), (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));

    // Handle the case where e.g. x is 1.99999988. x0 will be 1 and x1 will be 3 (1.99999988f + 1.0f = 3.0f).
    // If x1 - x0 > 1, x1--
    __m256 sub_mask = _mm256_cmp_ps(_mm256_sub_ps(x1, x0), _mm256_set1_ps(1.0f), _CMP_GT_OQ);
    x1 = _mm256_sub_ps(x1, _mm256_and_ps(_mm256_set1_ps(1.0f), sub_mask));
    sub_mask = _mm256_cmp_ps(_mm256_sub_ps(y1, y0), _mm256_set1_ps(1.0f), _CMP_GT_OQ);
    y1 = _mm256_sub_ps(y1, _mm256_and_ps(_mm256_set1_ps(1.0f), sub_mask));
    sub_mask = _mm256_cmp_ps(_mm256_sub_ps(z1, z0), _mm256_set1_ps(1.0f), _CMP_GT_OQ);
    z1 = _mm256_sub_ps(z1, _mm256_and_ps(_mm256_set1_ps(1.0f), sub_mask));

    // Smooth t.
    const __m256 dx0 = _mm256_sub_ps(x, x0);
    const __m256 dy0 = _mm256_sub_ps(y, y0);
    const __m256 dz0 = _mm256_sub_ps(z, z0);
    const __m256 tx = _mm256_mul_ps(dx0, _mm256_mul_ps(dx0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dx0, _mm256_set1_ps(2.0f)))));
    const __m256 ty = _mm256_mul_ps(dy0, _mm256_mul_ps(dy0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dy0, _mm256_set1_ps(2.0f)))));
    const __m256 tz = _mm256_mul_ps(dz0, _mm256_mul_ps(dz0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dz0, _mm256_set1_ps(2.0f)))));

    __m256 p0 = pnoise8_calc_gradient(x, y, z, x0, y0, z0);
    __m256 p1 = pnoise8_calc_gradient(x, y, z, x1, y0, z0);
    __m256 r0 = lerp8(p0, p1, tx);
    p0 = pnoise8_calc_gradient(x, y, z, x0, y1, z0);
    p1 = pnoise8_calc_gradient(x, y, z, x1, y1, z0);
    __m256 r1 = lerp8(p0, p1, tx);
    __m256 r2 = lerp8(r0, r1, ty);

    p0 = pnoise8_calc_gradient(x, y, z, x0, y0, z1);
    p1 = pnoise8_calc_gradient(x, y, z, x1, y0, z1);
    r0 = lerp8(p0, p1, tx);
    p0 = pnoise8_calc_gradient(x, y, z, x0, y1, z1);
    p1 = pnoise8_calc_gradient(x, y, z, x1, y1, z1);
    r1 = lerp8(p0, p1, tx);
    __m256 r3 = lerp8(r0, r1, ty);
    
    __m256 r4 = lerp8(r2, r3, tz);

    return r4;
}


void do_one_frame(struct MemoryArena* memory_arena)
{
    (void)memory_arena;

    static v3 cam_pos = {.m = {0.0f, 1.0f, 3.0f}};

#if 0


    static u8_m my_counter = 0xEE;
    my_counter += (u8_m)is_key_toggled_down(KB_N);
    u8 mc_idx = my_counter;

    u32_m num_grid_points = 8;
    v3* grid_points = (v3*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(v3) * 8);
    u8_m* grid_points_filled = (u8_m*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(u8) * 8);

    grid_points[0] = make_v3(0.0f, 0.0f, 0.0f);
    grid_points[1] = make_v3(2.0f, 0.0f, 0.0f);
    grid_points[2] = make_v3(0.0f, 0.0f, 2.0f);
    grid_points[3] = make_v3(2.0f, 0.0f, 2.0f);
    grid_points[4] = make_v3(0.0f, 2.0f, 0.0f);
    grid_points[5] = make_v3(2.0f, 2.0f, 0.0f);
    grid_points[6] = make_v3(0.0f, 2.0f, 2.0f);
    grid_points[7] = make_v3(2.0f, 2.0f, 2.0f);

    grid_points_filled[0] = mc_idx & (1 << 0);
    grid_points_filled[1] = mc_idx & (1 << 1);
    grid_points_filled[2] = mc_idx & (1 << 2);
    grid_points_filled[3] = mc_idx & (1 << 3);
    grid_points_filled[4] = mc_idx & (1 << 4);
    grid_points_filled[5] = mc_idx & (1 << 5);
    grid_points_filled[6] = mc_idx & (1 << 6);
    grid_points_filled[7] = mc_idx & (1 << 7);

    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));

    terrain_mesh->num_vertices = 24;
    for(u64_m i = 0; i < 3; i++)
    {
        terrain_mesh->vx[8*i + 0] = 0.0f;
        terrain_mesh->vx[8*i + 1] = 1.0f;
        terrain_mesh->vx[8*i + 2] = 2.0f;
        terrain_mesh->vx[8*i + 3] = 0.0f;
        terrain_mesh->vx[8*i + 4] = 2.0f;
        terrain_mesh->vx[8*i + 5] = 0.0f;
        terrain_mesh->vx[8*i + 6] = 1.0f;
        terrain_mesh->vx[8*i + 7] = 2.0f;

        terrain_mesh->vy[8*i + 0] = (f32)i;
        terrain_mesh->vy[8*i + 1] = (f32)i;
        terrain_mesh->vy[8*i + 2] = (f32)i;
        terrain_mesh->vy[8*i + 3] = (f32)i;
        terrain_mesh->vy[8*i + 4] = (f32)i;
        terrain_mesh->vy[8*i + 5] = (f32)i;
        terrain_mesh->vy[8*i + 6] = (f32)i;
        terrain_mesh->vy[8*i + 7] = (f32)i;

        terrain_mesh->vz[8*i + 0] = 0.0f;
        terrain_mesh->vz[8*i + 1] = 0.0f;
        terrain_mesh->vz[8*i + 2] = 0.0f;
        terrain_mesh->vz[8*i + 3] = 1.0f;
        terrain_mesh->vz[8*i + 4] = 1.0f;
        terrain_mesh->vz[8*i + 5] = 2.0f;
        terrain_mesh->vz[8*i + 6] = 2.0f;
        terrain_mesh->vz[8*i + 7] = 2.0f;
    }
    for(u64_m i = 0; i < terrain_mesh->num_vertices; i++)
    {
        u32 c = rand_u32(~(u32)(i));
        terrain_mesh->nx[i] = (f32)((c >> 16) & 0xFF) / 256.0f;
        terrain_mesh->ny[i] = (f32)((c >>  8) & 0xFF) / 256.0f;
        terrain_mesh->nz[i] = (f32)((c >>  0) & 0xFF) / 256.0f;
    }

    u32 num_tris = MARCHING_CUBES_NUM_TRIS[mc_idx];
    terrain_mesh->num_indices = num_tris * 3;
    for(u8_m i_tri = 0; i_tri < num_tris; i_tri++)
    {
        terrain_mesh->indices[3*i_tri + 0] = MARCHING_CUBES_TRIS_INDICES[mc_idx][i_tri][0];
        terrain_mesh->indices[3*i_tri + 1] = MARCHING_CUBES_TRIS_INDICES[mc_idx][i_tri][1];
        terrain_mesh->indices[3*i_tri + 2] = MARCHING_CUBES_TRIS_INDICES[mc_idx][i_tri][2];
    }

#elif 1

    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));

    s32 terrain_width = 64;
    s32 terrain_height = 64;
    s32 terrain_depth = 64;

    f32 x_freq = 0.06f;
    f32 y_freq = 0.06f;
    f32 z_freq = 0.06f;

    u32_m num_grid_points = 0;
    v3* grid_points = (v3*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(v3) * terrain_width * terrain_height * terrain_depth);
    u8_m* grid_points_filled = (u8_m*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(u8) * terrain_width * terrain_height * terrain_depth);

    u32_m out_num_indices = 0;
    u32_m* out_indices = terrain_mesh->indices;
    u32_m out_num_vertices = 0;
    f32_m* out_vx = terrain_mesh->vx;
    f32_m* out_vy = terrain_mesh->vy;
    f32_m* out_vz = terrain_mesh->vz;
    f32_m* out_nx = terrain_mesh->nx;
    f32_m* out_ny = terrain_mesh->ny;
    f32_m* out_nz = terrain_mesh->nz;


    // sample_x = round_neg_inf(cam_pos.x) - (f32)terrain_width * 0.5f + offset;
    const __m256 samples_x_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_round_ps(_mm256_set1_ps(cam_pos.x), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                _mm256_set1_ps((f32)terrain_width * 0.5f)
            ),
            _mm256_setr_ps(0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f)
        );
    // sample_y = round_neg_inf(cam_pos.y) - (f32)terrain_height * 0.5f + offset;
    const __m256 samples_y_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_round_ps(_mm256_set1_ps(cam_pos.y), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                _mm256_set1_ps((f32)terrain_height * 0.5f)
            ),
            _mm256_setr_ps(0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f)
        );
    // sample_z = round_neg_inf(cam_pos.z) - (f32)terrain_depth * 0.5f + offset;
    const __m256 samples_z_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_round_ps(_mm256_set1_ps(cam_pos.z), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                _mm256_set1_ps((f32)terrain_depth * 0.5f)
            ),
            _mm256_setr_ps(0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f)
        );
    for(s32_m i_y = 0; i_y < terrain_height; i_y++)
    {
        for(s32_m i_z = 0; i_z < terrain_depth; i_z++)
        {
            for(s32_m i_x = 0; i_x < terrain_width; i_x++)
            {
                const __m256 samples_x = _mm256_add_ps(_mm256_set1_ps((f32)i_x), samples_x_start);
                const __m256 samples_y = _mm256_add_ps(_mm256_set1_ps((f32)i_y), samples_y_start);
                const __m256 samples_z = _mm256_add_ps(_mm256_set1_ps((f32)i_z), samples_z_start);

                __m256 noise = pnoise8(
                        _mm256_mul_ps(samples_x, _mm256_set1_ps(x_freq)),
                        _mm256_mul_ps(samples_y, _mm256_set1_ps(y_freq)),
                        _mm256_mul_ps(samples_z, _mm256_set1_ps(z_freq))
                    );
                noise = _mm256_sub_ps(noise, _mm256_mul_ps(samples_y, _mm256_set1_ps(0.02f)));

                f32_m tris_x[5 * 3];
                f32_m tris_y[5 * 3];
                f32_m tris_z[5 * 3];
                u32 num_tris = marching_cube(tris_x, tris_y, tris_z, noise);

                ///////////////////////////////
                (void)num_grid_points;
                (void)grid_points;
                (void)grid_points_filled;
                //v3 sample_pos = {.x = (f32)sample_x, .y = (f32)sample_y, .z = (f32)sample_z };
                //grid_points[num_grid_points] = sample_pos;
                //grid_points_filled[num_grid_points] = mc_idx & 1;
                //num_grid_points++;
                ///////////////////////////////

                _Alignas(32) f32_m samples_x_arr[8];
                _Alignas(32) f32_m samples_y_arr[8];
                _Alignas(32) f32_m samples_z_arr[8];
                _mm256_store_ps(samples_x_arr, samples_x);
                _mm256_store_ps(samples_y_arr, samples_y);
                _mm256_store_ps(samples_z_arr, samples_z);

                for(u64_m i_tri = 0; i_tri < num_tris; i_tri++)
                {
                    v3 vert0 = { .x = tris_x[i_tri * 3 + 2], .y = tris_y[i_tri * 3 + 2], .z = tris_z[i_tri * 3 + 2] };
                    v3 vert1 = { .x = tris_x[i_tri * 3 + 1], .y = tris_y[i_tri * 3 + 1], .z = tris_z[i_tri * 3 + 1] };
                    v3 vert2 = { .x = tris_x[i_tri * 3 + 0], .y = tris_y[i_tri * 3 + 0], .z = tris_z[i_tri * 3 + 0] };

                    v3 normal = v3_normalize(v3_cross(v3_sub(vert2, vert0), v3_sub(vert1, vert0)));

                    v3 sample_offset = {.x = samples_x_arr[0], .y = samples_y_arr[0], .z = samples_z_arr[0] };
                    vert0 = v3_add(vert0, sample_offset);
                    vert1 = v3_add(vert1, sample_offset);
                    vert2 = v3_add(vert2, sample_offset);
                    
                    *out_vx++ = vert0.x;
                    *out_vy++ = vert0.y;
                    *out_vz++ = vert0.z;
                    *out_nx++ = normal.x;
                    *out_ny++ = normal.y;
                    *out_nz++ = normal.z;

                    *out_vx++ = vert1.x;
                    *out_vy++ = vert1.y;
                    *out_vz++ = vert1.z;
                    *out_nx++ = normal.x;
                    *out_ny++ = normal.y;
                    *out_nz++ = normal.z;

                    *out_vx++ = vert2.x;
                    *out_vy++ = vert2.y;
                    *out_vz++ = vert2.z;
                    *out_nx++ = normal.x;
                    *out_ny++ = normal.y;
                    *out_nz++ = normal.z;

                    out_num_vertices += 3;

                    *out_indices++ = out_num_indices++;
                    *out_indices++ = out_num_indices++;
                    *out_indices++ = out_num_indices++;
                }
            }
        }
    }
    terrain_mesh->num_vertices = out_num_vertices;
    terrain_mesh->num_indices = out_num_indices;

#else



    u32 terrain_width = 128;
    u32 terrain_vertex_stride = terrain_width + 1;
    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));
    terrain_mesh->num_vertices = (terrain_width + 1) * (terrain_width + 1);
    terrain_mesh->num_indices = 0;
    {
        f32 offset = 32.0f;
        f32 scale = 0.036f;
        __m256 accum_x = _mm256_fmadd_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(scale), _mm256_set1_ps(offset));
        __m256 accum_z = _mm256_set1_ps(offset);
        for(u64_m i_z = 0; i_z < terrain_vertex_stride; i_z++)
        {
            for(u64_m i_x = 0; i_x < terrain_vertex_stride - 8; i_x += 8)
            {
                __m256 v = pnoise8(accum_x, _mm256_set1_ps(0.0f), accum_z);

                f32_m storage8[8];
                _mm256_storeu_ps(storage8, v);

                for(u64_m lane = 0; lane < 8; lane++)
                {
                    f32 x = (f32)i_x - (f32)terrain_width * 0.5f;
                    f32 z = (f32)i_z - (f32)terrain_width * 0.5f;
                    f32 n = storage8[lane] * 0.5f + 0.5f;
                    terrain_mesh->vx[i_z * terrain_vertex_stride + i_x + lane] = x + (f32)lane;
                    terrain_mesh->vy[i_z * terrain_vertex_stride + i_x + lane] = n * 20.0f;
                    terrain_mesh->vz[i_z * terrain_vertex_stride + i_x + lane] = z;

                    u32 c = rand_u32((u32)(i_z * terrain_vertex_stride + i_x + lane));
                    terrain_mesh->nx[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >> 16) & 0xFF) / 256.0f;
                    terrain_mesh->ny[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >>  8) & 0xFF) / 256.0f;
                    terrain_mesh->nz[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >>  0) & 0xFF) / 256.0f;
                }

                accum_x = _mm256_add_ps(accum_x, _mm256_set1_ps(scale * 8.0f));
            }

            {
                __m256 v = pnoise8(accum_x, _mm256_set1_ps(0.0f), accum_z);
                f32_m storage8[8];
                _mm256_storeu_ps(storage8, v);
                for(u64_m i_x = terrain_vertex_stride & ~0b111, lane = 0; i_x < terrain_vertex_stride; i_x++, lane++)
                {
                    f32 x = (f32)i_x - (f32)terrain_width * 0.5f;
                    f32 z = (f32)i_z - (f32)terrain_width * 0.5f;
                    f32 n = storage8[lane] * 0.5f + 0.5f;
                    terrain_mesh->vx[i_z * terrain_vertex_stride + i_x] = x + (f32)lane;
                    terrain_mesh->vy[i_z * terrain_vertex_stride + i_x] = n * 20.0f;
                    terrain_mesh->vz[i_z * terrain_vertex_stride + i_x] = z;

                    terrain_mesh->nx[i_z * terrain_vertex_stride + i_x] = n;
                    terrain_mesh->ny[i_z * terrain_vertex_stride + i_x] = n;
                    terrain_mesh->nz[i_z * terrain_vertex_stride + i_x] = n;
                }
            }

            accum_z = _mm256_add_ps(accum_z, _mm256_set1_ps(scale));
            accum_x = _mm256_fmadd_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(scale), _mm256_set1_ps(offset));
        }

        for(u32_m i_z = 0; i_z < terrain_width; i_z++)
        {
            for(u32_m i_x = 0; i_x < terrain_width; i_x++)
            {
                u32 bl = (i_z + 1) * terrain_vertex_stride + i_x;
                u32 br = (i_z + 1) * terrain_vertex_stride + i_x + 1;
                u32 tl = (i_z + 0) * terrain_vertex_stride + i_x;
                u32 tr = (i_z + 0) * terrain_vertex_stride + i_x + 1;

                terrain_mesh->indices[terrain_mesh->num_indices++] = bl;
                terrain_mesh->indices[terrain_mesh->num_indices++] = br;
                terrain_mesh->indices[terrain_mesh->num_indices++] = tl;

                terrain_mesh->indices[terrain_mesh->num_indices++] = tl;
                terrain_mesh->indices[terrain_mesh->num_indices++] = br;
                terrain_mesh->indices[terrain_mesh->num_indices++] = tr;
            }
        }
    }
#endif

    {
        // Camera controls
        static f32_m cam_pitch_turns = 0.0f;
        cam_pitch_turns += (f32)is_key_down(KB_I) * 0.001f;
        cam_pitch_turns -= (f32)is_key_down(KB_K) * 0.001f;


        static f32_m cam_yaw_turns = 0.0f;
        cam_yaw_turns += (f32)is_key_down(KB_J) * 0.001f;
        cam_yaw_turns -= (f32)is_key_down(KB_L) * 0.001f;

        if(is_fps_mode())
        {
            s32_m mouse_screen_dx;
            s32_m mouse_screen_dy;
            get_mouse_delta(&mouse_screen_dx, &mouse_screen_dy);
            cam_pitch_turns -= (f32)(mouse_screen_dy) * 0.0005f;
            cam_yaw_turns -= (f32)(mouse_screen_dx) * 0.0005f;
        }


        mtx4x4 y_rot_mtx;
        make_y_axis_rotation_mtx(&y_rot_mtx, -cam_yaw_turns);

        mtx4x4 x_rot_mtx;
        make_x_axis_rotation_mtx(&x_rot_mtx, -cam_pitch_turns);

        mtx4x4 translation_mtx;
        make_translation_mtx(&translation_mtx, v3_scale(cam_pos, -1.0f));

        mtx4x4 world_to_cam_mtx_temp;
        mtx4x4_mul(&world_to_cam_mtx_temp, &y_rot_mtx, &translation_mtx);

        mtx4x4 world_to_cam_mtx;
        mtx4x4_mul(&world_to_cam_mtx, &x_rot_mtx, &world_to_cam_mtx_temp);

        v3 I_cam = {.m={world_to_cam_mtx.m[0*4 + 0], world_to_cam_mtx.m[0*4 + 1], world_to_cam_mtx.m[0*4 + 2]}};
        v3 J_cam = {.m={world_to_cam_mtx.m[1*4 + 0], world_to_cam_mtx.m[1*4 + 1], world_to_cam_mtx.m[1*4 + 2]}};
        v3 K_cam = {.m={world_to_cam_mtx.m[2*4 + 0], world_to_cam_mtx.m[2*4 + 1], world_to_cam_mtx.m[2*4 + 2]}};

        f32 speed = 1.0f;
        cam_pos = v3_add(cam_pos, v3_scale(K_cam, -speed * (f32)is_key_down(KB_W)));
        cam_pos = v3_add(cam_pos, v3_scale(K_cam,  speed * (f32)is_key_down(KB_S)));

        cam_pos = v3_add(cam_pos, v3_scale(I_cam,  speed * (f32)is_key_down(KB_D)));
        cam_pos = v3_add(cam_pos, v3_scale(I_cam, -speed * (f32)is_key_down(KB_A)));
        
        cam_pos = v3_add(cam_pos, v3_scale(J_cam,  speed * (f32)is_key_down(KB_SPACE)));
        cam_pos = v3_add(cam_pos, v3_scale(J_cam, -speed * (f32)is_key_down(KB_LCTRL)));

        /*
         * Derivation for 3D perspective projection matrix
         *        
         *                                                    /
         *                                                  /
         *                                                /
         *                                              /
         *                                            /
         *                                          /
         *                                        /
         *                                      /
         *                                    /
         *                                  /
         *                                /
         *                              /
         *                            /
         *                          /
         *                        /
         *                      / |               V
         *                    /   |             .>+---------------+
         *                  /     |       -----/ D|               |
         *         Y      /       |  ----/        |               |
         *         ^    /       --R-/             |               |
         *         |  /    ----/  |               |               |
         *         |/  ---/       |               |               |
         *  Z <----C../           |----> N        |               |
         *  F     / \             |               |               |
         *       /    \           |               +---------------+
         *      V       \         |
         *     X          \       |
         *                  \     |
         *                    \   |
         *                      \ |
         *                        \
         *                          \
         *                            \
         *                              \
         *                                \
         *                                  \
         *                                    \
         *                                      \
         *                                        \
         *                                          \
         *                                            \
         *                                              \
         *                                                \
         *                                                  \
         *                                                    \
         *        
         *        
         *        
         * C : 3D camera point (Assume the camera is at the origin - (0, 0)
         * F : Normalized camera forward vector (The camera looks along the -Z axis. For something to be seen, it must be more -Z than the near plane.)
         *     In camera space, this will be (0, 0, 1)
         * n : Camera's near plane distance. For something to be seen, it must have a Z coordinate < -n.
         * V : Vertex to be projected
         *
         * The goal is to intersect the ray from the origin to the vertex with the near plane.
         *
         * Q : Point along the ray (solve for intersection)
         * Q = C + unit(V)*t, but since C is just the origin
         * Q = unit(V)*t
         *
         * Find the plane equation:
         * N : The near plane normal (In camera space, this will be (0, 0, -1)
         * S : a point on the near plane.
         * S = N*n
         * 
         * New plane equation : (P - S) * N = 0
         * We want to find where a point on the ray is equal to 0, so plug in Q for P:
         * (Q - S) * N = 0
         *
         * Expand:
         *
         * (unit(V)*t - S) * N = 0
         *
         * Solve for t:
         *
         * unit(V)*t*N - S*N = 0
         * t = S*N / (unit(V)*N)
         *
         * Plug t back in to the ray equation:
         *
         * Q = unit(V)*t
         * Q = unit(V)*(S*N / (unit(V)*N))
         *
         * Find in terms of V
         *
         * Q = unit(V)*(S*N / (unit(V)*N))
         *
         * Q = unit(V)*S*N
         *     -----------
         *     (unit(V)*N)
         *
         * Q = unit(V)*(N*n)*N
         *     ---------------
         *       (unit(V)*N)
         *      
         * Q = N*n*N*unit(V)
         *     -------------
         *      (unit(V)*N)
         *
         * Q = N*n*N*(V / ||V*V||)
         *     -------------------
         *      ((V / ||V*V||)*N)
         *
         * Q = N*n*N*V
         *     -------
         *     (V*N)
         *
         * Q = N*N*n*V
         *     -------
         *     (V*N)
         *
         * In 3D:
         * Q = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
         *     -----------------------------------  *  V
         *       (V_x*N_x + V_y*N_y + V_z*N_z)
         *
         * Q_x = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
         *       -----------------------------------  *  V_x
         *         (V_x*N_x + V_y*N_y + V_z*N_z)
         *
         * Q_y = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
         *       -----------------------------------  *  V_y
         *         (V_x*N_x + V_y*N_y + V_z*N_z)
         *
         * Q_z = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
         *       -----------------------------------  *  V_z
         *         (V_x*N_x + V_y*N_y + V_z*N_z)
         * 
         * Assume our object has been translated to camera space. In this case, N = (0, 0, -1) (Right-handed coordinate system)
         * N_x = 0
         * N_y = 0
         * N_z = -1
         *
         * Q_x =   n
         *       ------ * V_x
         *       (-V_z)
         *
         * Q_y =   n
         *       ------ * V_y
         *       (-V_z)
         *
         * Q_z =   n
         *       ------ * V_z
         *       (-V_z)
         *
         *
         * So now we have the point Q in camera space where Q is V perspective projected onto the near plane.
         * Our goal is to find Q_p in NDC space. So, we need to divide X and Y by the camera width and height.
         *
         * C_w : camera width
         * C_h : camera height
         *
         * Q_px = Q_x / C_w = n / C_w
         *                    ------- * V_x
         *                    (-V_z)        
         *
         * Q_py = Q_y / C_h = n / C_h
         *                    ------- * V_y
         *                    (-V_z)
         *
         * Z should be between -1 and 1, so we need to divide by far plane - near plane. Use the vertex's Z coordinate instead of Q's Z coordinate
         * (Q is already projected and will have a constant Z, so we can't use that).
         * f : far plane dist
         *
         * Q_pz = (V_z - (-n)) * 2
         *        ----------------  -  1
         *           -f - (-n)
         *
         * Q_pz = (V_z + n) * 2
         *        -------------  -  1
         *           n - f
         *
         * Q_pz = V_z*2 + n*2
         *        -----------  -  1
         *           n - f
         *
         * Q_pz = V_z*2       n*2
         *        ------  +  -----  -  1
         *        n - f      n - f
         *
         * Q_pz =   2              n*2
         *        ----- * V_z  +  -----  -  1
         *        n - f           n - f
         *
         * https://www.desmos.com/calculator/frzetn7doc
         *       
         * Now, define as a matrix (keep in mind we will be dividing by the W component after matrix multiplication):
         *
         * | Q_px |   | n / C_w    0           0            0           |   |  V_x |   
         * | Q_py | = |   0      n / C_h       0            0           | * |  V_y |
         * | Q_pz |   |   0        0        2 / (n-f)    (n*2) / (n-f)  |   |  V_z |   
         * | Q_pw |   |   0        0          -1            0           |   | 1.0f |   
         *
         * | Q_px |   |         (n / C_w) * V_x         |
         * | Q_py | = |         (n / C_h) * V_y         |
         * | Q_pz |   | 2 / (n-f) * V_z + (n*2) / (n-f) |
         * | Q_pw |   |              -V_z               | <-- Will be dividing all the terms by -V_z
         * 
         * Dividing the depth by -V_z has the unfortunate consequence of reducing the NDC depth space.
         *
         */

        f32 aspect_ratio = get_screen_aspect_ratio();

        f32 n = 0.1f;
        f32 f = 1000.0f;
        f32 C_w = 0.125f;
        f32 C_h = C_w * aspect_ratio;
        mtx4x4 proj_mtx = {
            .m = {
                //    X         Y                Z                      W
                n / C_w,     0.0f,            0.0f,                  0.0f,
                   0.0f,  n / C_h,            0.0f,                  0.0f,
                   0.0f,     0.0f,  2.0f / (n - f),  (n * 2.0f) / (n - f),
                   0.0f,     0.0f,           -1.0f,                  0.0f,
            }
        };

        mtx4x4 camera_and_clip_mtx;
        mtx4x4_mul(&camera_and_clip_mtx, &proj_mtx, &world_to_cam_mtx);

        f32 pos[3] = {0};
        draw_Mesh_1M(terrain_mesh, &camera_and_clip_mtx, 1, pos + 0, pos + 1, pos + 2);

#if 0
        // Draw marching cube vertices
        {
            f32 circle_radius = 0.05f;
            v3 c_white = {.m={1.0f, 1.0f, 1.0f}};
            v3 c_black = {.m={0.0f, 0.0f, 0.0f}};

            for(u64_m i = 0; i < num_grid_points; i++)
            {
                debug_draw_sphere(&proj_mtx, &world_to_cam_mtx, grid_points[i], circle_radius, grid_points_filled[i] ? c_white : c_black);
            }
        }
#endif

        // Draw world basis.
        {
            v3 s = {.m={0.0f, 0.0f, 0.0f}};
            v3 e = {.m={1.0f, 0.0f, 0.0f}};
            debug_draw_line(&camera_and_clip_mtx, s, e, e);
            e.x = 0.0f;
            e.y = 1.0f;
            e.z = 0.0f;
            debug_draw_line(&camera_and_clip_mtx, s, e, e);
            e.x = 0.0f;
            e.y = 0.0f;
            e.z = 1.0f;
            debug_draw_line(&camera_and_clip_mtx, s, e, e);
        }
    }

}


