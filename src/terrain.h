
#pragma once

#include "common.h"
#include "math.h"

struct MemoryArena;

#define TERRAIN_MAX_CHUNKS 3000
#define CHUNK_DIM 16
#define CHUNK_POW 4
#define TERRAIN_CHUNK_MAX_VERTS 8192

INTERNAL inline s32 truncate_chunk(const s32 n)
{
    return n & ~15;
}

struct TerrainChunkGeometry
{
    u32 num_verts;
    f32 vx[TERRAIN_CHUNK_MAX_VERTS];
    f32 vy[TERRAIN_CHUNK_MAX_VERTS];
    f32 vz[TERRAIN_CHUNK_MAX_VERTS];
    f32 nx[TERRAIN_CHUNK_MAX_VERTS];
    f32 ny[TERRAIN_CHUNK_MAX_VERTS];
    f32 nz[TERRAIN_CHUNK_MAX_VERTS];

    _Static_assert(TERRAIN_CHUNK_MAX_VERTS <= u16_MAX, "TERRAIN_CHUNK_MAX_VERTS overflow.");
    u32 num_indices;
    u16 indices[TERRAIN_CHUNK_MAX_VERTS * 3];
};
struct Terrain
{
    u64 num_chunks;
    u64 chunks_key[TERRAIN_MAX_CHUNKS];
    u8 chunks_exponent[TERRAIN_MAX_CHUNKS];
    struct TerrainChunkGeometry chunks_geometry[TERRAIN_MAX_CHUNKS];
};

struct TerrainProgress
{
    s32 i;
};

INTERNAL u64 pack_voxel_key(const s32 x, const s32 y, const s32 z)
{
    ASSERT(abs_s32(x) < (1 << 19), "Voxel position will overflow when packing.");
    ASSERT(abs_s32(y) < (1 << 19), "Voxel position will overflow when packing.");
    ASSERT(abs_s32(z) < (1 << 19), "Voxel position will overflow when packing.");
    return ((const u64)z & 0xFFFFF) << 40ULL |
           ((const u64)y & 0xFFFFF) << 20ULL |
           ((const u64)x & 0xFFFFF);
}

INTERNAL void unpack_voxel_key(s32* x, s32* y, s32* z, const u64 k)
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
    _Alignas(32) s32 arr[8];
    _mm256_store_si256((__m256i*)arr, kv);
    *x = arr[1];
    *y = arr[3];
    *z = arr[5];
}

/*
INTERNAL inline struct TerrainChunk* get_chunk_packed(struct Terrain* terrain, const u64 key)
{
    u64 chunk_idx = (u64)-1;
    for(u64 i = 0; i < terrain->num_chunks; i++)
    {
        if(terrain->chunks_key[i] == key)
        {
            chunk_idx = i;
        }
    }
    return chunk_idx == (u64)-1 ? NULL : terrain->chunks + chunk_idx;
}

INTERNAL inline struct TerrainChunk* get_chunk(struct Terrain* terrain, const s32 x, const s32 y, const s32 z)
{
    const u64 key = pack_voxel_key(x, y, z);
    return get_chunk_packed(terrain, key);
}
*/

INTERNAL void reset_terrain(struct Terrain* terrain)
{
    terrain->num_chunks = 0;
}

INTERNAL void copy_terrain(struct Terrain* dst, const struct Terrain* src)
{
    dst->num_chunks = src->num_chunks;
    COPY(dst->chunks_key, src->chunks_key, src->num_chunks);
    COPY(dst->chunks_exponent, src->chunks_exponent, src->num_chunks);
    COPY(dst->chunks_geometry, src->chunks_geometry, src->num_chunks);
}

void generate_terrain(
        struct Terrain* terrain,
        struct TerrainProgress* progress,
        const f32 cam_pos_x,
        const f32 cam_pos_y,
        const f32 cam_pos_z,
        const u32 gen_max);
