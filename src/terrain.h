
#pragma once

#include "common.h"
#include "math.h"

struct MemoryArena;

#define TERRAIN_MAX_CHUNKS 3000
#define CHUNK_DIM 16
#define TERRAIN_CHUNK_MAX_VERTS 8192

INTERNAL inline s32 truncate_chunk(s32 n)
{
    return n & ~15;
}

struct TerrainChunk
{
    // Indexed by z * CHUNK_DIM*CHUNK_DIM + y*CHUNK_DIM + x
    _Alignas(32) f32_m m_val[CHUNK_DIM*CHUNK_DIM*CHUNK_DIM];
};
struct TerrainChunkGeometry
{
    u32_m num_verts;
    f32_m vx[TERRAIN_CHUNK_MAX_VERTS];
    f32_m vy[TERRAIN_CHUNK_MAX_VERTS];
    f32_m vz[TERRAIN_CHUNK_MAX_VERTS];
    f32_m nx[TERRAIN_CHUNK_MAX_VERTS];
    f32_m ny[TERRAIN_CHUNK_MAX_VERTS];
    f32_m nz[TERRAIN_CHUNK_MAX_VERTS];

    _Static_assert(TERRAIN_CHUNK_MAX_VERTS <= u16_MAX, "TERRAIN_CHUNK_MAX_VERTS overflow.");
    u32_m num_indices;
    u16_m indices[TERRAIN_CHUNK_MAX_VERTS * 3];
};
struct Terrain
{
    u64_m num_chunks;
    u64_m chunks_key[TERRAIN_MAX_CHUNKS];
    u8_m chunks_exponent[TERRAIN_MAX_CHUNKS];
    struct TerrainChunk chunks[TERRAIN_MAX_CHUNKS];
    struct TerrainChunkGeometry chunks_geometry[TERRAIN_MAX_CHUNKS];
};

struct TerrainProgress
{
    s32_m i;
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

INTERNAL inline struct TerrainChunk* get_chunk_packed(struct Terrain* terrain, u64 key)
{
    u64_m chunk_idx = (u64_m)-1;
    for(u64_m i = 0; i < terrain->num_chunks; i++)
    {
        if(terrain->chunks_key[i] == key)
        {
            chunk_idx = i;
        }
    }
    return chunk_idx == (u64_m)-1 ? NULL : terrain->chunks + chunk_idx;
}

INTERNAL inline struct TerrainChunk* get_chunk(struct Terrain* terrain, s32 x, s32 y, s32 z)
{
    u64 key = pack_voxel_key(x, y, z);
    return get_chunk_packed(terrain, key);
}

INTERNAL void reset_terrain(struct Terrain* terrain)
{
    terrain->num_chunks = 0;
}

INTERNAL void copy_terrain(struct Terrain* dst, const struct Terrain* src)
{
    dst->num_chunks = src->num_chunks;
    COPY(dst->chunks_key, src->chunks_key, src->num_chunks);
    COPY(dst->chunks_exponent, src->chunks_exponent, src->num_chunks);
    COPY(dst->chunks, src->chunks, src->num_chunks);
    COPY(dst->chunks_geometry, src->chunks_geometry, src->num_chunks);
}

void generate_terrain(
        struct Terrain* terrain,
        struct TerrainProgress* progress,
        f32 cam_pos_x,
        f32 cam_pos_y,
        f32 cam_pos_z,
        u32 gen_max);
