
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








// Thoughts on terrain generation:
//
// The goal here is to generate a flat-ish world with features that can form along all axis (that is to say, can cross back on itself in any axis).
//
// Ideas:
// * Generate a 2D heightmap terrain. Iterate over the generated chunks and apply a 3D modification to give overhands. 
//   This could work, but a key property we want to preserve is the ability to query any point in space and tell whether it's ground or air.
//   The 2D heightmap and 3D modification process would need to match the function used to query a point. The 2 different modes would be:
//   Generating terrain:
//   1. Iterate (x, z) points in a 2D grid and generate a y value. Add chunks.
//   2. Iterate chunks, for each point in the chunk, add or substract its noise values based on a secondary noise function. If the function bleeds into neighbors,
//      regenerate them. Repeat until all chunks are ground or air. Is this basically the surface walker idea?
//   Querying terrain:
//   1. Given (x, y, z), generate the height value from (x, z).
//   2. Add the secondary noise function given (x, y, z)
//   3. Check whether the noise function is <= 0.
//
// * Write a surface walking algorithm. Find a start point on the terrain and generate neighboring chunks containing the surface. Recursively generate terrain
//   outwards.
//   This would miss details like small floating islands. The Dijkstra-like "open-list" would get pretty massive.
//
// * Iterate over a (x, z) region. As you iterate, keep track of a "y-window" which has a variable min and max. As you iterate, the y-window will change and
//   tell you what possible chunks there are to generate. The y-window itself can be driven by a noies function.
//
//
// I think the y-window option is what I wanna do because it looks closest to brute force. If we go too far down the path of the surface walker and I decide I
// don't like, we're screwed.
//
// Are option 1 and 3 basically the same?
//
// Changing my mind. I want the terrain generation to be entirely driven by noise so it's simple to query. Forcing a range of Y values sounds hacky and restrictive.
// Let's try option 2 and see where it leads us. Then we can either pivot or simplify.
//


#if 0
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
#else
INTERNAL inline __m256 sample_terrain(const __m256 x, const __m256 y, const __m256 z)
{
    __m256 noise = pnoise8(
            _mm256_mul_ps(x, _mm256_set1_ps(0.03f)),
            _mm256_mul_ps(y, _mm256_set1_ps(0.03f)),
            _mm256_mul_ps(z, _mm256_set1_ps(0.03f)));
    //noise = _mm256_fmadd_ps(
    //        y,
    //        _mm256_set1_ps(0.03f),
    //        noise);
    noise = _mm256_fmadd_ps(
            y,
            _mm256_set1_ps(0.0075f),
            noise);
    return noise;
}

INTERNAL inline f32 sample_terrain_scalar(f32 s_x, f32 s_y, f32 s_z)
{
    _Alignas(32) f32_m arr[8];
    _mm256_store_ps(
            arr,
            sample_terrain(
                _mm256_set1_ps(s_x),
                _mm256_set1_ps(s_y),
                _mm256_set1_ps(s_z)));
    return arr[0];
}
#endif

u8 generate_chunk(struct TerrainChunk* chunk, s32 chunk_x, s32 chunk_y, s32 chunk_z)
{
    u8_m res_mask = 0;
    u8_m is_first_iter = 1;

    __m256 x8 = _mm256_setr_ps(
        (f32)(chunk_x + 0),
        (f32)(chunk_x + 1),
        (f32)(chunk_x + 2),
        (f32)(chunk_x + 3),
        (f32)(chunk_x + 4),
        (f32)(chunk_x + 5),
        (f32)(chunk_x + 6),
        (f32)(chunk_x + 7)
    );
    __m256 y8 = _mm256_set1_ps((f32)chunk_y);
    __m256 z8 = _mm256_set1_ps((f32)chunk_z);

    for(s32_m i_y = 0; i_y < CHUNK_DIM; i_y++)
    {
        for(s32_m i_z = 0; i_z < CHUNK_DIM; i_z++)
        {
            for(s32_m i_x = 0; i_x < CHUNK_DIM; i_x += 8)
            {
                _Alignas(32) const __m256 noise = sample_terrain(x8, y8, z8);
                _mm256_store_ps(chunk->m_val + i_z * CHUNK_DIM*CHUNK_DIM + i_y*CHUNK_DIM + i_x, noise);

                u8 mask = (u8)_mm256_movemask_ps(noise);

                res_mask = res_mask ^ mask;
                res_mask = is_first_iter ? 0 : res_mask;
                is_first_iter = 0;

                x8 = _mm256_add_ps(x8, _mm256_set1_ps(8.0f));
            }

            x8 = _mm256_setr_ps(
                (f32)(chunk_x + 0),
                (f32)(chunk_x + 1),
                (f32)(chunk_x + 2),
                (f32)(chunk_x + 3),
                (f32)(chunk_x + 4),
                (f32)(chunk_x + 5),
                (f32)(chunk_x + 6),
                (f32)(chunk_x + 7)
            );
            z8 = _mm256_add_ps(z8, _mm256_set1_ps(1.0f));
        }

        z8 = _mm256_set1_ps((f32)chunk_z);
        y8 = _mm256_add_ps(y8, _mm256_set1_ps(1.0f));
    }

    return res_mask > 0;
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

    u32_m num_chunks = 0;
    for(s32_m chunk_y = -64; chunk_y < 64; chunk_y += CHUNK_DIM)
    {
        for(s32_m chunk_z = -64; chunk_z < 64; chunk_z += CHUNK_DIM)
        {
            for(s32_m chunk_x = -64; chunk_x < 64; chunk_x += CHUNK_DIM)
            {
                ASSERT(num_chunks < TERRAIN_MAX_CHUNKS, "Overflow max chunks.");

                terrain->chunks_key[num_chunks] = pack_voxel_key(chunk_x, chunk_y, chunk_z);
                terrain->chunks_exponent[num_chunks] = 1;
                struct TerrainChunk* chunk = terrain->chunks + num_chunks;
                u8 add_chunk = generate_chunk(chunk, chunk_x, chunk_y, chunk_z);

                num_chunks += add_chunk;
            }
        }
    }
    terrain->num_chunks = num_chunks;



    for(u64_m i_chunk = 0; i_chunk < terrain->num_chunks; i_chunk++)
    {
        struct TerrainChunkGeometry* chunk_geo = terrain->chunks_geometry + i_chunk;
        f32_m* dst_vx = chunk_geo->vx;
        f32_m* dst_vy = chunk_geo->vy;
        f32_m* dst_vz = chunk_geo->vz;
        f32_m* dst_nx = chunk_geo->nx;
        f32_m* dst_ny = chunk_geo->ny;
        f32_m* dst_nz = chunk_geo->nz;
        u16_m* dst_indices = chunk_geo->indices;
        u32_m num_verts = 0;
        u16_m num_indices = 0;

        struct TerrainChunk* chunk_000 = terrain->chunks + i_chunk;
        s32_m chunk_x;
        s32_m chunk_y;
        s32_m chunk_z;
        unpack_voxel_key(&chunk_x, &chunk_y, &chunk_z, terrain->chunks_key[i_chunk]);

        struct TerrainChunk* chunk_001 = get_chunk(terrain, chunk_x + 1, chunk_y + 0, chunk_z + 0);
        struct TerrainChunk* chunk_010 = get_chunk(terrain, chunk_x + 0, chunk_y + 1, chunk_z + 0);
        struct TerrainChunk* chunk_011 = get_chunk(terrain, chunk_x + 1, chunk_y + 1, chunk_z + 0);
        struct TerrainChunk* chunk_100 = get_chunk(terrain, chunk_x + 0, chunk_y + 0, chunk_z + 1);
        struct TerrainChunk* chunk_101 = get_chunk(terrain, chunk_x + 1, chunk_y + 0, chunk_z + 1);
        struct TerrainChunk* chunk_110 = get_chunk(terrain, chunk_x + 0, chunk_y + 1, chunk_z + 1);
        struct TerrainChunk* chunk_111 = get_chunk(terrain, chunk_x + 1, chunk_y + 1, chunk_z + 1);

        // TODO(mfritz) neighbors
        for(s32_m i_z0 = 0; i_z0 < CHUNK_DIM - 1; i_z0++)
        {
            for(s32_m i_y0 = 0; i_y0 < CHUNK_DIM - 1; i_y0++)
            {
                for(s32_m i_x0 = 0; i_x0 < CHUNK_DIM - 1; i_x0++)
                {
                    f32 x0 = (f32)chunk_x + i_x0;
                    f32 y0 = (f32)chunk_y + i_y0;
                    f32 z0 = (f32)chunk_z + i_z0;
                    f32 x1 = (f32)chunk_x + i_x0 + 1;
                    f32 y1 = (f32)chunk_y + i_y0 + 1;
                    f32 z1 = (f32)chunk_z + i_z0 + 1;

                    u8 x_on_edge = i_x0 + 1 == CHUNK_DIM;
                    u8 y_on_edge = i_y0 + 1 == CHUNK_DIM;
                    u8 z_on_edge = i_z0 + 1 == CHUNK_DIM;
                    s32 i_x1 = x_on_edge ? 0 : i_x0 + 1;
                    s32 i_y1 = y_on_edge ? 0 : i_y0 + 1;
                    s32 i_z1 = z_on_edge ? 0 : i_z0 + 1;

                    struct TerrainChunk* chunk_table[] = {
                        chunk_000,
                        chunk_001,
                        chunk_010,
                        chunk_011,

                        chunk_100,
                        chunk_101,
                        chunk_110,
                        chunk_111,
                    };

                    struct TerrainChunk* chunk_for_sample_000 = chunk_table[0 * 4 +         0 * 2 +         0];
                    struct TerrainChunk* chunk_for_sample_001 = chunk_table[0 * 4 +         0 * 2 + x_on_edge];
                    struct TerrainChunk* chunk_for_sample_010 = chunk_table[0 * 4 + y_on_edge * 2 +         0];
                    struct TerrainChunk* chunk_for_sample_011 = chunk_table[0 * 4 + y_on_edge * 2 + x_on_edge];

                    struct TerrainChunk* chunk_for_sample_100 = chunk_table[z_on_edge * 4 +         0 * 2 +         0];
                    struct TerrainChunk* chunk_for_sample_101 = chunk_table[z_on_edge * 4 +         0 * 2 + x_on_edge];
                    struct TerrainChunk* chunk_for_sample_110 = chunk_table[z_on_edge * 4 + y_on_edge * 2 +         0];
                    struct TerrainChunk* chunk_for_sample_111 = chunk_table[z_on_edge * 4 + y_on_edge * 2 + x_on_edge];

                    _Alignas(32) f32_m vals[8] = {};
                    // _mm256_store_ps(vals, sample_terrain(
                    //     _mm256_setr_ps(x0, x1, x1, x0, x0, x1, x1, x0),
                    //     _mm256_setr_ps(y0, y0, y0, y0, y1, y1, y1, y1),
                    //     _mm256_setr_ps(z0, z0, z1, z1, z0, z0, z1, z1)
                    // ));
                    (void)x1;
                    (void)y1;
                    (void)z1;

                    vals[0] = 
                        chunk_for_sample_000
                        ? chunk_for_sample_000->m_val[i_z0 * CHUNK_DIM*CHUNK_DIM + i_y0*CHUNK_DIM + i_x0]
                        : vals[0];
                    vals[1] =
                        chunk_for_sample_001
                        ? chunk_for_sample_001->m_val[i_z0 * CHUNK_DIM*CHUNK_DIM + i_y0*CHUNK_DIM + i_x1]
                        : vals[1];
                    vals[2] =
                        chunk_for_sample_101
                        ? chunk_for_sample_101->m_val[i_z1 * CHUNK_DIM*CHUNK_DIM + i_y0*CHUNK_DIM + i_x1]
                        : vals[2];
                    vals[3] =
                        chunk_for_sample_100
                        ? chunk_for_sample_100->m_val[i_z1 * CHUNK_DIM*CHUNK_DIM + i_y0*CHUNK_DIM + i_x0]
                        : vals[3];

                    vals[4] =
                        chunk_for_sample_010
                        ? chunk_for_sample_010->m_val[i_z0 * CHUNK_DIM*CHUNK_DIM + i_y1*CHUNK_DIM + i_x0]
                        : vals[4];
                    vals[5] =
                        chunk_for_sample_011
                        ? chunk_for_sample_011->m_val[i_z0 * CHUNK_DIM*CHUNK_DIM + i_y1*CHUNK_DIM + i_x1]
                        : vals[5];
                    vals[6] =
                        chunk_for_sample_111
                        ? chunk_for_sample_111->m_val[i_z1 * CHUNK_DIM*CHUNK_DIM + i_y1*CHUNK_DIM + i_x1]
                        : vals[6];
                    vals[7] =
                        chunk_for_sample_110
                        ? chunk_for_sample_110->m_val[i_z1 * CHUNK_DIM*CHUNK_DIM + i_y1*CHUNK_DIM + i_x0]
                        : vals[7];

                    f32_m tris_x[5 * 3];
                    f32_m tris_y[5 * 3];
                    f32_m tris_z[5 * 3];
                    u32 num_tris = marching_cube(tris_x, tris_y, tris_z, _mm256_load_ps(vals));

                    for(u64_m i_tri = 0; i_tri < num_tris; i_tri++)
                    {
                        v3 v_0 = make_v3(tris_x[i_tri * 3 + 0], tris_y[i_tri * 3 + 0], tris_z[i_tri * 3 + 0]);
                        v3 v_1 = make_v3(tris_x[i_tri * 3 + 1], tris_y[i_tri * 3 + 1], tris_z[i_tri * 3 + 1]);
                        v3 v_2 = make_v3(tris_x[i_tri * 3 + 2], tris_y[i_tri * 3 + 2], tris_z[i_tri * 3 + 2]);

                        v3 n = v3_normalize(v3_cross(v3_sub(v_1, v_0), v3_sub(v_2, v_0)));

                        *dst_vx++ = tris_x[i_tri * 3 + 0] + x0;
                        *dst_vy++ = tris_y[i_tri * 3 + 0] + y0;
                        *dst_vz++ = tris_z[i_tri * 3 + 0] + z0;
                        *dst_nx++ = n.x;
                        *dst_ny++ = n.y;
                        *dst_nz++ = n.z;
                        *dst_indices++ = num_indices;
                        num_indices++;
                        ASSERT(num_indices < u16_MAX, "Vertex index overflow.");
                        num_verts++;

                        *dst_vx++ = tris_x[i_tri * 3 + 1] + x0;
                        *dst_vy++ = tris_y[i_tri * 3 + 1] + y0;
                        *dst_vz++ = tris_z[i_tri * 3 + 1] + z0;
                        *dst_nx++ = n.x;
                        *dst_ny++ = n.y;
                        *dst_nz++ = n.z;
                        *dst_indices++ = num_indices;
                        num_indices++;
                        ASSERT(num_indices < u16_MAX, "Vertex index overflow.");
                        num_verts++;

                        *dst_vx++ = tris_x[i_tri * 3 + 2] + x0;
                        *dst_vy++ = tris_y[i_tri * 3 + 2] + y0;
                        *dst_vz++ = tris_z[i_tri * 3 + 2] + z0;
                        *dst_nx++ = n.x;
                        *dst_ny++ = n.y;
                        *dst_nz++ = n.z;
                        *dst_indices++ = num_indices;
                        num_indices++;
                        ASSERT(num_indices < u16_MAX, "Vertex index overflow.");
                        num_verts++;

                        ASSERT(num_verts < TERRAIN_CHUNK_MAX_VERTS, "Exceeded max verts per chunk.");
                        ASSERT(num_indices < TERRAIN_CHUNK_MAX_VERTS * 3, "Exceeded max indices per chunk.");
                    }
                }
            }
        }

        chunk_geo->num_verts = num_verts;
        chunk_geo->num_indices = num_indices;
    }

}




