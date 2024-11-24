
#include "terrain.h"
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

INTERNAL inline __m256 sample_terrain(const __m256 x, const __m256 y, const __m256 z)
{
    __m256 layer0 = pnoise8(
            _mm256_mul_ps(x, _mm256_set1_ps(0.02f)),
            _mm256_mul_ps(y, _mm256_set1_ps(0.02f)),
            _mm256_mul_ps(z, _mm256_set1_ps(0.02f))
        );
    //layer0 = _mm256_sub_ps(layer0, _mm256_mul_ps(y, _mm256_set1_ps(0.02f)));
    __m256 layer1 = pnoise8(
            _mm256_mul_ps(x, _mm256_set1_ps(0.01f)),
            _mm256_mul_ps(y, _mm256_set1_ps(0.01f)),
            _mm256_mul_ps(z, _mm256_set1_ps(0.01f))
        );
    layer1 = _mm256_sub_ps(layer1, _mm256_mul_ps(y, _mm256_set1_ps(0.01f)));
    __m256 layer2 = pnoise8(
            _mm256_mul_ps(x, _mm256_set1_ps(0.0015f)),
            _mm256_mul_ps(y, _mm256_set1_ps(0.0025f)),
            _mm256_mul_ps(z, _mm256_set1_ps(0.0015f))
        );
    layer2 = _mm256_sub_ps(layer2, _mm256_mul_ps(y, _mm256_set1_ps(0.0f)));

    __m256 result = _mm256_fmadd_ps(
        _mm256_set1_ps(0.1f),
        layer0,
        _mm256_fmadd_ps(
            _mm256_set1_ps(0.2f),
            layer1,
            _mm256_mul_ps(
                _mm256_set1_ps(2.0f),
                layer2
            )
        )
    );
    return result;
}

INTERNAL inline struct TerrainChunk* get_chunk(struct Terrain* terrain, s32 x, s32 y, s32 z)
{
    u64 key = pack_voxel_key(x, y, z);
    u64_m chunk_idx = (u64_m)-1;
    for(u64_m i = 0; i < terrain->num_chunks; i++)
    {
        if(terrain->chunk_key[i] == key)
        {
            chunk_idx = i;
        }
    }
    return chunk_idx == (u64_m)-1 ? NULL : terrain->chunks + chunk_idx;
}

void generate_terrain(
        struct Terrain* terrain,
        struct TerrainProgress* progress,
        f32 cam_pos_x,
        f32 cam_pos_y,
        f32 cam_pos_z,
        u32 gen_max)
{
    (void)terrain;
    (void)progress;
    (void)cam_pos_x;
    (void)cam_pos_y;
    (void)cam_pos_z;
    (void)gen_max;

    terrain->num_chunks = 0;

    __m256 s_x8 = _mm256_setr_ps(
            -64.0f + 0.0f,
            -64.0f + 1.0f,
            -64.0f + 2.0f,
            -64.0f + 3.0f,
            -64.0f + 4.0f,
            -64.0f + 5.0f,
            -64.0f + 6.0f,
            -64.0f + 7.0f);
    __m256 s_z8 = _mm256_setr_ps(
            -64.0f,
            -64.0f,
            -64.0f,
            -64.0f,
            -64.0f,
            -64.0f,
            -64.0f,
            -64.0f);
    for(s32_m s_z = -128; s_z < 128; s_z++)
    {
        for(s32_m s_x = -128; s_x < 128; s_x += 8)
        {
            __m256 noise = pnoise8(
                    _mm256_mul_ps(s_x8, _mm256_set1_ps(0.01f)),
                    _mm256_set1_ps(5.0f),
                    _mm256_mul_ps(s_z8, _mm256_set1_ps(0.01f)));
            noise = _mm256_mul_ps(noise, _mm256_set1_ps(70.0f));
            _Alignas(32) f32_m noise_arr[8];
            _mm256_store_ps(noise_arr, noise);

            for(s32_m i_arr = 0; i_arr < 8; i_arr++)
            {
                s32 chunk_x = truncate_chunk(s_x + i_arr);
                s32 chunk_y = truncate_chunk((s32)noise_arr[i_arr]);
                s32 chunk_z = truncate_chunk(s_z);
                struct TerrainChunk* chunk = get_chunk(terrain, chunk_x, chunk_y, chunk_z);
                if(!chunk)
                {
                    u64 num_chunks = terrain->num_chunks;
                    ASSERT(num_chunks < TERRAIN_MAX_CHUNKS, "Chunk overflow");
                    terrain->chunk_key[num_chunks] = pack_voxel_key(chunk_x, chunk_y, chunk_z);
                    terrain->chunk_exponent[num_chunks] = 1;
                    terrain->num_chunks++;
                    chunk = terrain->chunks + num_chunks;
                }

                s32 x = s_x + i_arr - chunk_x;
                s32 y = (s32)noise_arr[i_arr] - chunk_y;
                s32 z = s_z - chunk_z;

                for(s32_m i_y = 0; i_y < y; i_y++)
                {
                    chunk->m_filled[z * CHUNK_DIM*CHUNK_DIM + i_y*CHUNK_DIM + x] = 1;
                }
                for(s32_m i_y = y; i_y < CHUNK_DIM; i_y++)
                {
                    chunk->m_filled[z * CHUNK_DIM*CHUNK_DIM + i_y*CHUNK_DIM + x] = 0;
                }
            }

            s_x8 = _mm256_add_ps(s_x8, _mm256_set1_ps(8.0f));
        }

        s_x8 = _mm256_setr_ps(
            -64.0f + 0.0f,
            -64.0f + 1.0f,
            -64.0f + 2.0f,
            -64.0f + 3.0f,
            -64.0f + 4.0f,
            -64.0f + 5.0f,
            -64.0f + 6.0f,
            -64.0f + 7.0f);
        s_z8 = _mm256_add_ps(s_z8, _mm256_set1_ps(1.0f));
    }

}




