
#include "common.h"
#include "platform.h"
#include "graphics.h"
#include "math.h"


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





// https://paulbourke.net/geometry/polygonise/
typedef struct {
   v3 p[3];
} TRIANGLE;

typedef struct {
   v3 p[8];
   f32_m val[8];
} GRIDCELL;
/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
v3 VertexInterp(f32 isolevel, v3 p1, v3 p2, f32 valp1, f32 valp2)
{
    if(absf(isolevel-valp1) < 0.00001f)
    {
        return p1;
    }
    if(absf(isolevel-valp2) < 0.00001f)
    {
        return p2;
    }
    if(absf(valp1-valp2) < 0.00001f)
    {
        return p1;
    }
    f32 mu = (isolevel - valp1) / (valp2 - valp1);
    v3 p;
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);
    return p;
}
int Polygonise(GRIDCELL grid, f32 isolevel, TRIANGLE *triangles)
{
    u32 edgeTable[256]={
        0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
        0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
        0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
        0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
        0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
        0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
        0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
        0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
        0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
        0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
        0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
        0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
        0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
        0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
        0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
        0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
        0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
        0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
        0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
        0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
        0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
    s32 triTable[256][16] =
    {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

    s32_m cubeindex = 0;
    if (grid.val[0] < isolevel) cubeindex |= 1;
    if (grid.val[1] < isolevel) cubeindex |= 2;
    if (grid.val[2] < isolevel) cubeindex |= 4;
    if (grid.val[3] < isolevel) cubeindex |= 8;
    if (grid.val[4] < isolevel) cubeindex |= 16;
    if (grid.val[5] < isolevel) cubeindex |= 32;
    if (grid.val[6] < isolevel) cubeindex |= 64;
    if (grid.val[7] < isolevel) cubeindex |= 128;

    /* Cube is entirely in/out of the surface */
    if (edgeTable[cubeindex] == 0)
        return(0);

    /* Find the vertices where the surface intersects the cube */
    v3 vertlist[12];
    if (edgeTable[cubeindex] & 1)    vertlist[0]  = VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
    if (edgeTable[cubeindex] & 2)    vertlist[1]  = VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
    if (edgeTable[cubeindex] & 4)    vertlist[2]  = VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
    if (edgeTable[cubeindex] & 8)    vertlist[3]  = VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
    if (edgeTable[cubeindex] & 16)   vertlist[4]  = VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
    if (edgeTable[cubeindex] & 32)   vertlist[5]  = VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
    if (edgeTable[cubeindex] & 64)   vertlist[6]  = VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
    if (edgeTable[cubeindex] & 128)  vertlist[7]  = VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
    if (edgeTable[cubeindex] & 256)  vertlist[8]  = VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
    if (edgeTable[cubeindex] & 512)  vertlist[9]  = VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
    if (edgeTable[cubeindex] & 1024) vertlist[10] = VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
    if (edgeTable[cubeindex] & 2048) vertlist[11] = VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

    /* Create the triangle */
    int num_tris = 0;
    for (int i=0; triTable[cubeindex][i] != -1; i += 3)
    {
        triangles[num_tris].p[0] = vertlist[triTable[cubeindex][i  ]];
        triangles[num_tris].p[1] = vertlist[triTable[cubeindex][i+1]];
        triangles[num_tris].p[2] = vertlist[triTable[cubeindex][i+2]];
        num_tris++;
    }

    return num_tris;
}


void do_one_frame(struct MemoryArena* memory_arena)
{
    (void)memory_arena;


    /*
     * 1. Create a table mapping 8-bit int -> triangulation
     * 2. Create a test grid with perlin noise values
     * 3. Display
     */

#if 1


    GRIDCELL grid = {0};
    grid.p[0].x = 0.0f;  grid.p[0].y = 0.0f;  grid.p[0].z = 0.0f;
    grid.p[1].x = 1.0f;  grid.p[1].y = 0.0f;  grid.p[1].z = 0.0f;
    grid.p[2].x = 1.0f;  grid.p[2].y = 0.0f;  grid.p[2].z = 1.0f;
    grid.p[3].x = 0.0f;  grid.p[3].y = 0.0f;  grid.p[3].z = 1.0f;
    grid.p[4].x = 0.0f;  grid.p[4].y = 1.0f;  grid.p[4].z = 0.0f;
    grid.p[5].x = 1.0f;  grid.p[5].y = 1.0f;  grid.p[5].z = 0.0f;
    grid.p[6].x = 1.0f;  grid.p[6].y = 1.0f;  grid.p[6].z = 1.0f;
    grid.p[7].x = 0.0f;  grid.p[7].y = 1.0f;  grid.p[7].z = 1.0f;
    for(u32_m iter = 0; iter < 256; iter++)
    {
        for(u32_m i = 0; i < 8; i++)
        {
            grid.val[i] = ((1 << i) & iter) ? 1.0f : -1.0f;
        }

        TRIANGLE out_tris[512];
        int num_tris = Polygonise(grid, 0.0f, out_tris);

        // Dump tris

        //char buf[1024];
        //snprintf(buf, sizeof(buf), "(%f, %f, %f)\n", out_tris[0].p[0].x, out_tris[0].p[0].y, out_tris[0].p[0].z);

        (void)num_tris;
        (void)out_tris;
        int b = 0;
        b++;
        (void)b;
    }

    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));
    (void)terrain_mesh;
    /*
    terrain_mesh->num_vertices = num_tris * 3;
    terrain_mesh->num_indices = num_tris * 3;

    for(s32_m i_tri = 0; i_tri < num_tris; i_tri++)
    {
        terrain_mesh->vx[3*i_tri + 0] = out_tris[i_tri].p[0].x;
        terrain_mesh->vy[3*i_tri + 0] = out_tris[i_tri].p[0].y;
        terrain_mesh->vz[3*i_tri + 0] = out_tris[i_tri].p[0].z;

        terrain_mesh->vx[3*i_tri + 1] = out_tris[i_tri].p[1].x;
        terrain_mesh->vy[3*i_tri + 1] = out_tris[i_tri].p[1].y;
        terrain_mesh->vz[3*i_tri + 1] = out_tris[i_tri].p[1].z;

        terrain_mesh->vx[3*i_tri + 2] = out_tris[i_tri].p[2].x;
        terrain_mesh->vy[3*i_tri + 2] = out_tris[i_tri].p[2].y;
        terrain_mesh->vz[3*i_tri + 2] = out_tris[i_tri].p[2].z;

        u32 c = rand_u32(~(u32)(i_tri));
        terrain_mesh->nx[3*i_tri + 0] = (f32)((c >> 16) & 0xFF) / 256.0f;
        terrain_mesh->ny[3*i_tri + 0] = (f32)((c >>  8) & 0xFF) / 256.0f;
        terrain_mesh->nz[3*i_tri + 0] = (f32)((c >>  0) & 0xFF) / 256.0f;

        terrain_mesh->nx[3*i_tri + 1] = (f32)((c >> 16) & 0xFF) / 256.0f;
        terrain_mesh->ny[3*i_tri + 1] = (f32)((c >>  8) & 0xFF) / 256.0f;
        terrain_mesh->nz[3*i_tri + 1] = (f32)((c >>  0) & 0xFF) / 256.0f;

        terrain_mesh->nx[3*i_tri + 2] = (f32)((c >> 16) & 0xFF) / 256.0f;
        terrain_mesh->ny[3*i_tri + 2] = (f32)((c >>  8) & 0xFF) / 256.0f;
        terrain_mesh->nz[3*i_tri + 2] = (f32)((c >>  0) & 0xFF) / 256.0f;

        terrain_mesh->indices[3*i_tri + 0] = 3*i_tri + 0;
        terrain_mesh->indices[3*i_tri + 1] = 3*i_tri + 1;
        terrain_mesh->indices[3*i_tri + 2] = 3*i_tri + 2;
    }
    */

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

        static v3 cam_pos = {.m = {0.0f, 1.0f, 3.0f}};

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

        f32 speed = 0.5f;
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


