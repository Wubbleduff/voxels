
#pragma once

#include <immintrin.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Math

INTERNAL inline f32 abs_f32(f32 a)
{
    return a < 0.0f ? -a : a;
} 

INTERNAL inline s32 abs_s32(s32 a)
{
    return a < 0 ? -a : a;
} 

INTERNAL inline f32 round_neg_inf(f32 a)
{
    return _mm_cvtss_f32(_mm_round_ps(_mm_set1_ps(a), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
}

#define H_PI 1.57079637f
#define PI 3.14159274f
#define TAU 6.28318548f
#define INV_PI 0.318309873f
#define INV_H_PI 0.636620f
#define INV_TAU 0.159154937f

typedef struct mtx4x4Tag
{
    union
    {
        _Alignas(32) f32_m m[16];
        __m128 v[4];
    };
} mtx4x4;

typedef struct v2Tag
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
        struct
        {
            f32_m x;
            f32_m y;
        };
    };
} v2;

typedef struct v3Tag
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
        struct
        {
            f32_m x;
            f32_m y;
            f32_m z;
        };
    };
} v3;

typedef struct v4Tag
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
        struct
        {
            f32_m x;
            f32_m y;
            f32_m z;
            f32_m w;
        };
    };
} v4;

INTERNAL inline v3 make_v3(f32 x, f32 y, f32 z)
{
    v3 r;
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}

INTERNAL inline v3 v3_zero()
{
    v3 r;
    r.x = 0.0f;
    r.y = 0.0f;
    r.z = 0.0f;
    return r;
}

INTERNAL inline f32 v3_dot(v3 a, v3 b)
{
    // | a.x | a.y | a.z | - |
    // | b.x | b.y | b.z | - |

    v3 r;

    //     00        01        10
    // | a.x*b.x | a.y*b.y | a.z*b.z | - |
    __m128 p = _mm_mul_ps(a.v, b.v);

    // | a.x*b.x | a.y*b.y | a.z*b.z | - |
    // +
    // | a.y*b.y | a.z*b.z | a.x*b.x | - |
    r.v = _mm_add_ps(p, _mm_shuffle_ps(p, p, 0b11001001));

    // | a.x*b.x + a.y*b.y | a.y*b.y + a.z*b.z | a.z*b.z + a.x*b.x | - |
    // +
    // |           a.z*b.z |           a.x*b.x |           a.y*b.y | - |
    r.v = _mm_add_ps(r.v, _mm_shuffle_ps(p, p, 0b11010010));

    return r.m[0];
}

INTERNAL inline v3 v3_cross(v3 a, v3 b)
{
    v3 r;
    r.v = _mm_fmsub_ps(
        _mm_shuffle_ps(a.v, a.v, 0b11001001),
        _mm_shuffle_ps(b.v, b.v, 0b11010010),
        _mm_mul_ps(
            _mm_shuffle_ps(a.v, a.v, 0b11010010),
            _mm_shuffle_ps(b.v, b.v, 0b11001001)
        )
    );
    return r;
}

INTERNAL inline v3 v3_add(v3 a, v3 b)
{
    v3 r;
    r.v = _mm_add_ps(a.v, b.v);
    return r;
}

INTERNAL inline v3 v3_sub(v3 a, v3 b)
{
    v3 r;
    r.v = _mm_sub_ps(a.v, b.v);
    return r;
}

INTERNAL inline v3 v3_scale(v3 a, f32 b)
{
    v3 r;
    r.v = _mm_mul_ps(a.v, _mm_set1_ps(b));
    return r;
}

INTERNAL inline v3 v3_normalize(v3 a)
{
    // TODO(mfritz) No need to go to scalar here.
    v3 r = a;
    f32 len_sqr = v3_dot(a, a);
    f32 len = _mm_cvtss_f32(_mm_sqrt_ss(_mm_set1_ps(len_sqr)));
    r = v3_scale(r, 1.0f / len);
    return r;
}

static inline __m256 lerp8(__m256 a, __m256 b, __m256 t)
{
    return _mm256_add_ps(
            _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1.0f), t), a),
            _mm256_mul_ps(t, b)
            );
}

INTERNAL inline __m256 approx_sin8(__m256 x)
{
    const __m256 pi       = _mm256_set1_ps(PI);
    const __m256 h_pi     = _mm256_set1_ps(H_PI);
    const __m256 inv_pi   = _mm256_set1_ps(INV_PI);
    const __m256 inv_h_pi = _mm256_set1_ps(INV_H_PI);
    const __m256 inv_tau = _mm256_set1_ps(INV_TAU);
    
    // Range reduce to [0, 2*pi] to check if we need to negate result.
    __m256 sign;
    {
        __m256 xx = _mm256_mul_ps(x, inv_tau);
        xx = _mm256_sub_ps(xx, _mm256_round_ps(xx, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
        __m256 mask = _mm256_cmp_ps(xx, _mm256_set1_ps(0.5f), _CMP_NLT_UQ);
        sign = _mm256_and_ps(mask, _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));
    }
    
    // From now on x is positive
    x = _mm256_and_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)));
    
    // Range reduce to [0, pi] to check if we need to flip X.
    __m256 flip_mask;
    {
        __m256 xx = _mm256_mul_ps(x, inv_pi);
        xx = _mm256_sub_ps(xx, _mm256_round_ps(xx, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
        flip_mask = _mm256_cmp_ps(xx, _mm256_set1_ps(0.5f), _CMP_NLT_UQ);
    }
    
    // Range reduce to [0, pi/2] for evaluation
    x = _mm256_mul_ps(x, inv_h_pi);
    x = _mm256_sub_ps(x, _mm256_round_ps(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    x = _mm256_mul_ps(x,  h_pi);
    
    __m256 flip_x = _mm256_fmadd_ps(x, _mm256_set1_ps(-1.0f), h_pi);
    x = _mm256_blendv_ps(x, flip_x, flip_mask);
    
    // x = x - x^3 / 6
    __m256 eval0 =  _mm256_fnmadd_ps(_mm256_mul_ps(_mm256_mul_ps(x, x), x), _mm256_set1_ps(1.0f / 6.0f), x);
    
    // x = 1/384*(pi - 2*x)^4 - 1/8*(pi - 2x)^2 + 1
    __m256 t = _mm256_fnmadd_ps(x, _mm256_set1_ps(2.0f), pi);
    __m256 t2 = _mm256_mul_ps(t, t);
    __m256 t4 = _mm256_mul_ps(t2, t2);
    __m256 eval1 = _mm256_mul_ps(_mm256_set1_ps(1.0f/384.0f), t4);
    eval1 = _mm256_fnmadd_ps(_mm256_set1_ps(1.0f/8.0f), t2, eval1);
    eval1 = _mm256_add_ps(eval1, _mm256_set1_ps(1.0f));
    
    __m256 eval_blend_mask = _mm256_cmp_ps(x, _mm256_set1_ps(0.6403f), _CMP_LT_OQ);
    __m256 eval = _mm256_blendv_ps(eval1, eval0, eval_blend_mask);
    
    eval = _mm256_or_ps(eval, sign);
    
    return eval;
}

INTERNAL inline __m256 approx_cos8(__m256 x)
{
    const __m256 h_pi = _mm256_set1_ps(0.5f * 3.14159265f);
    return approx_sin8(_mm256_sub_ps(x, h_pi));
}

// 4x4 matrix multiply : r = a * b.
// NOTE: Matrices are assumed to be row-major.
INTERNAL inline void mtx4x4_mul(mtx4x4* r, mtx4x4* a, mtx4x4* b)
{
    /*
     *
     * | a00  a01  a02  a03 |   | b00  b01  b02  b03 |   
     * | a10  a11  a12  a13 | * | b10  b11  b12  b13 | = 
     * | a20  a21  a22  a23 |   | b20  b21  b22  b23 |   
     * | a30  a31  a32  a33 |   | b30  b31  b32  b33 |   
     *
     * | a00*b00 + a01*b10 + a02*b20 + a03*b30    a00*b01 + a01*b11 + a02*b21 + a03*b31    a00*b02 + a01*b12 + a02*b22 + a03*b32    a00*b03 + a01*b13 + a02*b23 + a03*b33 |
     * | a10*b00 + a11*b10 + a12*b20 + a13*b30    a10*b01 + a11*b11 + a12*b21 + a13*b31    a10*b02 + a11*b12 + a12*b22 + a13*b32    a10*b03 + a11*b13 + a12*b23 + a13*b33 |
     * | a20*b00 + a21*b10 + a22*b20 + a23*b30    a20*b01 + a21*b11 + a22*b21 + a23*b31    a20*b02 + a21*b12 + a22*b22 + a23*b32    a20*b03 + a21*b13 + a22*b23 + a23*b33 |
     * | a30*b00 + a31*b10 + a32*b20 + a33*b30    a30*b01 + a31*b11 + a32*b21 + a33*b31    a30*b02 + a31*b12 + a32*b22 + a33*b32    a30*b03 + a31*b13 + a32*b23 + a33*b33 |
     *
     */

    /*
    Reference
    for(u64_m i = 0; i < 4; i++)
    {
        for(u64_m j = 0; j < 4; j++)
        {
            r[i*4 + j] =
                a[i*4 + 0] * b[0*4 + j] + 
                a[i*4 + 1] * b[1*4 + j] + 
                a[i*4 + 2] * b[2*4 + j] + 
                a[i*4 + 3] * b[3*4 + j];
        }
    }
    */

    for(u64_m row = 0; row < 4; row++)
    {
        const __m128 v = _mm_fmadd_ps(
            _mm_broadcast_ss(a->m + row*4 + 0),
            _mm_loadu_ps(b->m + 0*4 + 0),
            _mm_fmadd_ps(
                _mm_broadcast_ss(a->m + row*4 + 1),
                _mm_loadu_ps(b->m + 1*4 + 0),
                _mm_fmadd_ps(
                    _mm_broadcast_ss(a->m + row*4 + 2),
                    _mm_loadu_ps(b->m + 2*4 + 0),
                    _mm_mul_ps(
                        _mm_broadcast_ss(a->m + row*4 + 3),
                        _mm_loadu_ps(b->m + 3*4 + 0)
                    )
                )
            )
        );
        _mm_storeu_ps(r->m + row*4, v);
    }
}

INTERNAL inline void make_x_axis_rotation_mtx(mtx4x4* r, f32 turns)
{
    __m256 a_v = approx_sin8(_mm256_setr_ps(TAU * turns, TAU * turns + H_PI, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
    _Alignas(32) f32_m a[8];
    _mm256_store_ps(a, a_v);
    f32 sin_a = a[0];
    f32 cos_a = a[1];

    r->m[0] = 1.0f;   r->m[1] = 0.0f;   r->m[2] = 0.0f;    r->m[3] = 0.0f;
    r->m[4] = 0.0f;   r->m[5] = cos_a;  r->m[6] = -sin_a;  r->m[7] = 0.0f;
    r->m[8] = 0.0f;   r->m[9] = sin_a;  r->m[10] = cos_a;  r->m[11] = 0.0f;
    r->m[12] = 0.0f;  r->m[13] = 0.0f;  r->m[14] = 0.0f;   r->m[15] = 1.0f;
}

INTERNAL inline void make_y_axis_rotation_mtx(mtx4x4* r, f32 turns)
{
    __m256 a_v = approx_sin8(_mm256_setr_ps(TAU * turns, TAU * turns + H_PI, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
    _Alignas(32) f32_m a[8];
    _mm256_store_ps(a, a_v);
    f32 sin_a = a[0];
    f32 cos_a = a[1];

    r->m[0] = cos_a;   r->m[1] = 0.0f;   r->m[2] = sin_a;   r->m[3] = 0.0f;
    r->m[4] = 0.0f;    r->m[5] = 1.0f;   r->m[6] = 0.0f;    r->m[7] = 0.0f;
    r->m[8] = -sin_a;  r->m[9] = 0.0f;   r->m[10] = cos_a;  r->m[11] = 0.0f;
    r->m[12] = 0.0f;   r->m[13] = 0.0f;  r->m[14] = 0.0f;   r->m[15] = 1.0f;
}

INTERNAL inline void make_translation_mtx(mtx4x4* r, v3 v)
{
    r->m[0]  = 1.0f;   r->m[1] = 0.0f;    r->m[2] = 0.0f;   r->m[3] = v.m[0];
    r->m[4]  = 0.0f;   r->m[5] = 1.0f;    r->m[6] = 0.0f;   r->m[7] = v.m[1];
    r->m[8]  = 0.0f;   r->m[9] = 0.0f;   r->m[10] = 1.0f;  r->m[11] = v.m[2];
    r->m[12] = 0.0f;  r->m[13] = 0.0f;   r->m[14] = 0.0f;  r->m[15] = 1.0f;
}

INTERNAL inline u32 rand_u32(u32_m n)
{
    n ^= n << 13;
    n ^= n >> 17;
    n ^= n << 5;
    return n;
}

// NOTE: Should not initialize this with the result of rand_u32(). That will cause the 8 lanes to produce the same random numbers in sequence repeatedly.
INTERNAL inline __m256i rand8_u32(__m256i n)
{
    n = _mm256_xor_si256(n, _mm256_slli_epi32(n, 13));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(182376581));
    n = _mm256_xor_si256(n, _mm256_srli_epi32(n, 17));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(783456103));
    n = _mm256_xor_si256(n, _mm256_slli_epi32(n, 5));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(53523));
    return n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
