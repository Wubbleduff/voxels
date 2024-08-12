
#pragma once

#include "common.h"
#include "math.h"

struct MemoryArena;

#define TERRAIN_MAX_NUM_LOD 6
#define TERRAIN_MAX_NUM_VOXELS 65536
#define LOD_REGION_DIM 128

// Sparse storage for terrain voxels
struct TerrainVoxels
{
    u64_m num_voxels;
    u64_m keys[TERRAIN_MAX_NUM_VOXELS];

    u64_m num_tris[TERRAIN_MAX_NUM_VOXELS];
    f32_m tris_x[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
    f32_m tris_y[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
    f32_m tris_z[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
    f32_m normals_x[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
    f32_m normals_y[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
    f32_m normals_z[TERRAIN_MAX_NUM_VOXELS * 5 * 3];
};
struct Terrain
{
    struct TerrainVoxels terrain_lod[TERRAIN_MAX_NUM_LOD];
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
    for(u64_m i = 0; i < ARRAY_COUNT(terrain->terrain_lod); i++)
    {
        terrain->terrain_lod[i].num_voxels = 0;
    }
}

INTERNAL void copy_terrain(struct Terrain* dst, const struct Terrain* src)
{
    for(u64_m i = 0; i < ARRAY_COUNT(src->terrain_lod); i++)
    {
        u64 num_voxels = src->terrain_lod[i].num_voxels;
        dst->terrain_lod[i].num_voxels = num_voxels;
        ARRAY_COPY(dst->terrain_lod[i].keys, src->terrain_lod[i].keys, num_voxels);
        ARRAY_COPY(dst->terrain_lod[i].num_tris, src->terrain_lod[i].num_tris, num_voxels);
        ARRAY_COPY(dst->terrain_lod[i].tris_x, src->terrain_lod[i].tris_x, num_voxels * 5 * 3);
        ARRAY_COPY(dst->terrain_lod[i].tris_y, src->terrain_lod[i].tris_y, num_voxels * 5 * 3);
        ARRAY_COPY(dst->terrain_lod[i].tris_z, src->terrain_lod[i].tris_z, num_voxels * 5 * 3);
        ARRAY_COPY(dst->terrain_lod[i].normals_x, src->terrain_lod[i].normals_x, num_voxels * 5 * 3);
        ARRAY_COPY(dst->terrain_lod[i].normals_y, src->terrain_lod[i].normals_y, num_voxels * 5 * 3);
        ARRAY_COPY(dst->terrain_lod[i].normals_z, src->terrain_lod[i].normals_z, num_voxels * 5 * 3);
    }
}

void generate_terrain(
        struct Terrain* terrain,
        f32 cam_pos_x,
        f32 cam_pos_y,
        f32 cam_pos_z);
