
#pragma once

#include "common.h"
#include "math.h"

struct MemoryArena;

#define TERRAIN_MAX_CHUNKS 65536
#define CHUNK_DIM 16

INTERNAL inline s32 truncate_chunk(s32 n)
{
    return n & ~15;
}

struct TerrainChunk
{
    // 16*16*16 bitarray
    //u8 m_filled[512];
    // Indexed by z * CHUNK_DIM*CHUNK_DIM + y*CHUNK_DIM + x
    u8_m m_filled[16*16*16];
};
struct Terrain
{
    u64_m num_chunks;
    u64_m chunk_key[TERRAIN_MAX_CHUNKS];
    u8_m chunk_exponent[TERRAIN_MAX_CHUNKS];
    struct TerrainChunk chunks[TERRAIN_MAX_CHUNKS];
};

struct TerrainProgress
{
    int b;
};

INTERNAL u64 pack_voxel_key(s32 x, s32 y, s32 z)
{
    ASSERT(abs_s32(x) < (1 << 19), "Voxel position will overflow when packing.");
    ASSERT(abs_s32(y) < (1 << 19), "Voxel position will overflow when packing.");
    ASSERT(abs_s32(z) < (1 << 19), "Voxel position will overflow when packing.");
    return ((u64)z & 0xFFFFF) << 40ULL |
           ((u64)y & 0xFFFFF) << 20ULL |
           ((u64)x & 0xFFFFF);
}

INTERNAL void unpack_voxel_key(s32_m* x, s32_m* y, s32_m* z, u64 k)
{
    __m256i kv = _mm256_set1_epi64x(k);
    kv = _mm256_sllv_epi64(
        kv,
        _mm256_setr_epi64x(44, 24, 4, 0)
    );
    kv = _mm256_srav_epi32(
        kv,
        _mm256_setr_epi32(0, 12, 0, 12, 0, 12, 0, 0)
    );
    _Alignas(32) s32_m arr[8];
    _mm256_store_si256((__m256i*)arr, kv);
    *x = arr[1];
    *y = arr[3];
    *z = arr[5];
}

INTERNAL void reset_terrain(struct Terrain* terrain)
{
    (void)terrain;
}

INTERNAL void copy_terrain(struct Terrain* dst, const struct Terrain* src)
{
    (void)dst;
    (void)src;
}

void generate_terrain(
        struct Terrain* terrain,
        struct TerrainProgress* progress,
        f32 cam_pos_x,
        f32 cam_pos_y,
        f32 cam_pos_z,
        u32 gen_max);
