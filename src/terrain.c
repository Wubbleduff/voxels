
#include "terrain.h"
#include "math.h"
#include "marching_cubes.h"

#include <immintrin.h>

// Goal:
// Generate terrain with view distance 262,144 m
//
// Need a generated region of 524288^3 = 1.44115188075855872 * 10^17
// To simplify terrain, will reduce LOD as distance increase from player.
//
// Example:
// LOD0_WIDTH = 8
// MAX_LOD = 3
//
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *       *       *       *       *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                      *       *       *       *       *       *       *       *       *                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *   *   *   *   *   *   *   *   *       *       *               *               *     
//                                                                      |                                                                     
//                                                      *   *   *   *   *   *   *   *   *                                                     
//                                                                      |                                                                     
//                                      *       *       *   *   * * * * * * * * *   *   *       *       *                                     
//                                                              * * * * * * * * *                                                             
//                                                      *   *   * * * * * * * * *   *   *                                                     
//                                                              * * * * * * * * *                                                             
//  ----*---------------*---------------*-------*-------*---*---*-*-*-*-*-*-*-*-*---*---*-------*-------*---------------*---------------*-------------------
//                                                              * * * * * * * * *                                                             
//                                                      *   *   * * * * * * * * *   *   *                                                     
//                                                              * * * * * * * * *                                                             
//                                      *       *       *   *   * * * * * * * * *   *   *       *       *                                     
//                                                                      |                                                                     
//                                                      *   *   *   *   *   *   *   *   *                                                     
//                                                                      |                                                                     
//      *               *               *       *       *   *   *   *   *   *   *   *   *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                      *       *       *       *       *       *       *       *       *                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *       *       *       *       *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                    9 8 7 6 5 4 3 2 1 | 1 2 3 4 5 6 7 8 9                                                   
// 
// In 3D:
// Each LOD is 2x wider / sparser than the last
// (just store all samples in each layer for simplicity)
// # samples = 8^3 + 8^3 + 8^3 + 8^3 = 2048
//
// Example:
// LOD0_WIDTH = 64
// MAX_LOD = 18
// # samples = 18 * 64^3 = 4718592 (~5 million)
//

#define LOD0_WIDTH 128
#define MAX_LOD 18

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

// Input:
// * LOD 0 width (static u32)
// * Max LOD (static u32)
// * Camera pos (3x s32)
// Output:
// * List of tris (vertices, normals, indices)
// * TODO Acceleration spatial lookup structure for finding terrain
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
        u8 lod)
{
    u32_m num_vertices = 0;
    u32_m num_indices = 0;

    f32 x_freq = 0.02f;
    f32 y_freq = 0.02f;
    f32 z_freq = 0.02f;
    f32 y_bias = 0.02f;

    f32 region_width = LOD0_WIDTH;
    f32 region_scale = (float)(1 << (u32)lod);


    // sample_x = round_neg_inf(cam_pos.x / region_scale) - (f32)region_width * 0.5f * region_scale + offset;
    const __m256 samples_x_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_mul_ps(
                    _mm256_round_ps(_mm256_set1_ps((float)cam_pos_x / region_scale), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                    _mm256_set1_ps(region_scale)
                ),
                _mm256_set1_ps((f32)region_width * 0.5f * region_scale)
            ),
            _mm256_setr_ps(0.0f, region_scale, region_scale, 0.0f, 0.0f, region_scale, region_scale, 0.0f)
        );
    // sample_y = round_neg_inf(cam_pos.y / region_scale) - (f32)region_width * 0.5f * region_scale + offset;
    const __m256 samples_y_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_mul_ps(
                    _mm256_round_ps(_mm256_set1_ps((float)cam_pos_y / region_scale), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                    _mm256_set1_ps(region_scale)
                ),
                _mm256_set1_ps((f32)region_width * 0.5f * region_scale)
            ),
            _mm256_setr_ps(0.0f, 0.0f, 0.0f, 0.0f, region_scale, region_scale, region_scale, region_scale)
        );
    // sample_z = round_neg_inf(cam_pos.z / region_scale) - (f32)region_width * 0.5f * region_scale + offset;
    const __m256 samples_z_start = _mm256_add_ps(
            _mm256_sub_ps(
                _mm256_mul_ps(
                    _mm256_round_ps(_mm256_set1_ps((float)cam_pos_z / region_scale), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC),
                    _mm256_set1_ps(region_scale)
                ),
                _mm256_set1_ps((f32)region_width * 0.5f * region_scale)
            ),
            _mm256_setr_ps(0.0f, 0.0f, region_scale, region_scale, 0.0f, 0.0f, region_scale, region_scale)
        );
    for(s32_m i_y = 0; i_y < region_width; i_y++)
    {
        for(s32_m i_z = 0; i_z < region_width; i_z++)
        {
            for(s32_m i_x = 0; i_x < region_width; i_x++)
            {
                const __m256 samples_x = _mm256_add_ps(_mm256_set1_ps((f32)i_x * region_scale), samples_x_start);
                const __m256 samples_y = _mm256_add_ps(_mm256_set1_ps((f32)i_y * region_scale), samples_y_start);
                const __m256 samples_z = _mm256_add_ps(_mm256_set1_ps((f32)i_z * region_scale), samples_z_start);

                __m256 noise = pnoise8(
                        _mm256_mul_ps(samples_x, _mm256_set1_ps(x_freq)),
                        _mm256_mul_ps(samples_y, _mm256_set1_ps(y_freq)),
                        _mm256_mul_ps(samples_z, _mm256_set1_ps(z_freq))
                    );
                noise = _mm256_sub_ps(noise, _mm256_mul_ps(samples_y, _mm256_set1_ps(y_bias)));

                f32_m tris_x[5 * 3];
                f32_m tris_y[5 * 3];
                f32_m tris_z[5 * 3];
                u32 num_tris = marching_cube(tris_x, tris_y, tris_z, noise);

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

                    vert0 = v3_scale(vert0, region_scale);
                    vert1 = v3_scale(vert1, region_scale);
                    vert2 = v3_scale(vert2, region_scale);

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

                    num_vertices += 3;

                    *out_indices++ = num_indices++;
                    *out_indices++ = num_indices++;
                    *out_indices++ = num_indices++;
                }
            }
        }
    }

    *out_num_vertices = num_vertices;
    *out_num_indices = num_indices;
}




