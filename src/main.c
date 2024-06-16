
#include "common.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <immintrin.h>

#include <GL/GL.h>
#include "wglext.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Globals and defines

HANDLE g_log_file;
struct OpenGLState* g_opengl_state;
struct InputState* g_input_state;
s64_m g_clock_freq;

// TODO(mfritz) Should not be global
struct DrawData* g_draw_data;

#define ENGINE_FRAME_DURATION_US 16666LL / 2LL

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

typedef struct
{
    // Last byte should always be 0 so we can trivially cast to C-string.
    _Alignas(32) char buf[32];
} StringBuf32;

static inline const char* StringBuf32_to_cstr(const StringBuf32* s)
{
    return s->buf;
}

static inline StringBuf32 u32_to_StringBuf32(u32_m n)
{
    StringBuf32 result;
    memset(result.buf, 0, sizeof(result));
    u64_m size = 0;
    while(n && size < 31)
    {
        u32 digit = n % 10;
        result.buf[size++] = (char)('0' + digit);
        n /= 10;
    }
    for(u64_m i = 0; i < size / 2; i++)
    {
        char tmp = result.buf[i];
        result.buf[i] = result.buf[size - i - 1];
        result.buf[size - i - 1] = tmp;
    }
    return result;
}


static inline u64 cstr_len(const char* s)
{
    u64_m n = 0;
    while(*s++) n++;
    return n;
}

static inline void assert_fn(const char* file, int line, u64 c, const char* msg)
{
    if(!c)
    {
        const StringBuf32 line_str = u32_to_StringBuf32(line);
        const char* line_cstr = StringBuf32_to_cstr(&line_str);
        WriteFile(g_log_file, file, (int)cstr_len(file), NULL, NULL);
        WriteFile(g_log_file, ":", 1, NULL, NULL);
        WriteFile(g_log_file, line_cstr, (int)cstr_len(line_cstr), NULL, NULL);
        WriteFile(g_log_file, "  ", 1, NULL, NULL);
        WriteFile(g_log_file, msg, (int)cstr_len(msg), NULL, NULL);
        FlushFileBuffers(g_log_file);
        DebugBreak();
    }
}
#define ASSERT(c, msg) assert_fn(__FILE__, __LINE__, (u64)(c), (msg))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Memory

struct MemoryArena
{
    u8* const end;
    u8* const start;
    u8_m* cur;
};

MAYBE_UNUSED static inline struct MemoryArena memory_arena_init(void* base, u64 cap)
{
    ASSERT(base, "Memory arena init with NULL base.");
    struct MemoryArena result = {
        .start = base,
        .end = (u8*)base + cap,
        .cur = (u8_m*)base,
    };
    return result;
}

MAYBE_UNUSED static inline void memory_arena_reset(struct MemoryArena* arena)
{
    arena->cur = (u8_m*)arena->start;
}

MAYBE_UNUSED static inline void* memory_arena_allocate(const char* file, s32 line, struct MemoryArena* arena, u64 req_bytes)
{
    (void)file;
    (void)line;
    // Align to 64 byte boundary.
    u64 bytes = (req_bytes + 63ULL) & ~63ULL;
    ASSERT(arena->cur + bytes <= arena->end, "Memory arena allocation overflow.");
    void* result = arena->cur;
    arena->cur += bytes;
    return result;
}

MAYBE_UNUSED static inline void* memory_arena_allocate_zeroed(const char* file, s32 line, struct MemoryArena* arena, u64 req_bytes)
{
    (void)file;
    (void)line;
    void* result = memory_arena_allocate(file, line, arena, req_bytes);
    memset(result, 0, req_bytes);
    return result;
}

#define MEMORY_ARENA_ALLOCATE(arena, req_bytes) memory_arena_allocate(__FILE__, __LINE__, (arena), (req_bytes))
#define MEMORY_ARENA_ALLOCATE_ZEROED(arena, req_bytes) memory_arena_allocate_zeroed(__FILE__, __LINE__, (arena), (req_bytes))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Math

#define H_PI 1.57079637f
#define PI 3.14159274f
#define TAU 6.28318548f
#define INV_PI 0.318309873f
#define INV_H_PI 0.636620f
#define INV_TAU 0.159154937f

typedef struct
{
    union
    {
        _Alignas(32) f32_m m[16];
        __m128 v[4];
    };
} mtx4x4_m;
typedef const mtx4x4_m mtx4x4;

typedef struct
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
} v2_m;
typedef const v2_m v2;

typedef struct
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
} v3_m;
typedef const v3_m v3;

typedef struct
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
} v4_m;
typedef const v4_m v4;

MAYBE_UNUSED static inline f32 v3_dot(v3 a, v3 b)
{
    // | a.x | a.y | a.z | - |
    // | b.x | b.y | b.z | - |

    v3_m r;

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

MAYBE_UNUSED static inline v3 v3_cross(v3 a, v3 b)
{
    v3_m r;
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

MAYBE_UNUSED static inline v3 v3_add(v3 a, v3 b)
{
    v3_m r;
    r.v = _mm_add_ps(a.v, b.v);
    return r;
}

MAYBE_UNUSED static inline v3 v3_scale(v3 a, f32 b)
{
    v3_m r;
    r.v = _mm_mul_ps(a.v, _mm_set1_ps(b));
    return r;
}

MAYBE_UNUSED static inline v3 v3_normalize(v3 a)
{
    // TODO(mfritz) No need to go to scalar here.
    v3_m r = a;
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

MAYBE_UNUSED static inline __m256 approx_sin8(__m256 x)
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

MAYBE_UNUSED static inline __m256 approx_cos8(__m256 x)
{
    const __m256 h_pi = _mm256_set1_ps(0.5f * 3.14159265f);
    return approx_sin8(_mm256_sub_ps(x, h_pi));
}

// 4x4 matrix multiply : r = a * b.
// NOTE: Matrices are assumed to be row-major.
MAYBE_UNUSED static inline void mtx4x4_mul(mtx4x4_m* r, mtx4x4* a, mtx4x4* b)
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

MAYBE_UNUSED static inline void make_x_axis_rotation_mtx(mtx4x4_m* r, f32 turns)
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

MAYBE_UNUSED static inline void make_y_axis_rotation_mtx(mtx4x4_m* r, float turns)
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

MAYBE_UNUSED static inline void make_translation_mtx(mtx4x4_m* r, v3 v)
{
    r->m[0]  = 1.0f;   r->m[1] = 0.0f;    r->m[2] = 0.0f;   r->m[3] = v.m[0];
    r->m[4]  = 0.0f;   r->m[5] = 1.0f;    r->m[6] = 0.0f;   r->m[7] = v.m[1];
    r->m[8]  = 0.0f;   r->m[9] = 0.0f;   r->m[10] = 1.0f;  r->m[11] = v.m[2];
    r->m[12] = 0.0f;  r->m[13] = 0.0f;   r->m[14] = 0.0f;  r->m[15] = 1.0f;
}

MAYBE_UNUSED static inline u32 rand_u32(u32_m n)
{
    n ^= n << 13;
    n ^= n >> 17;
    n ^= n << 5;
    return n;
}

// NOTE: Should not initialize this with the result of rand_u32(). That will cause the 8 lanes to produce the same random numbers in sequence repeatedly.
MAYBE_UNUSED static inline __m256i rand8_u32(__m256i n)
{
    n = _mm256_xor_si256(n, _mm256_slli_epi32(n, 13));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(182376581));
    n = _mm256_xor_si256(n, _mm256_srli_epi32(n, 17));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(783456103));
    n = _mm256_xor_si256(n, _mm256_slli_epi32(n, 5));
    n = _mm256_mullo_epi32(n, _mm256_set1_epi32(53523));
    return n;
}


MAYBE_UNUSED static inline __m256 pnoise8_calc_gradient(__m256 x, __m256 y, __m256 z, __m256 vx, __m256 vy, __m256 vz)
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

MAYBE_UNUSED static inline __m256 pnoise8(const __m256 x, const __m256 y, const __m256 z)
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL / WGL

// Definitions copied from glad library.
typedef signed long long int GLsizeiptr;
typedef signed long long int GLintptr;
typedef char GLchar;
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_ARRAY_BUFFER 0x8892
#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#define GL_CLAMP_TO_EDGE 0x812F

// https://www.khronos.org/opengl/wiki/Load_OpenGL_Functions
void *load_gl_fn(HMODULE opengl32_dll_module, const char *name)
{
    void* p = (void*)wglGetProcAddress(name);
    if(p == 0 || (p == (void*)0x1) || (p == (void*)0x2) || (p == (void*)0x3) || (p == (void*)-1) )
    {
        p = (void *)GetProcAddress(opengl32_dll_module, name);
    }

    ASSERT(p != NULL, "Could not load OpenGL function.");
    ASSERT(p != (void*)1, "Could not load OpenGL function.");
    ASSERT(p != (void*)2, "Could not load OpenGL function.");
    ASSERT(p != (void*)-1, "Could not load OpenGL function.");

    return p;
}

typedef HGLRC (*fnptr_wglCreateContextAttribsARB)(HDC hDC, HGLRC hshareContext, const int *attribList);
typedef void (*fnptr_glGenVertexArrays)(GLsizei n, GLuint *arrays);
typedef GLenum (*fnptr_glGetError)();
typedef void (*fnptr_glBindVertexArray)(GLuint array);
typedef void (*fnptr_glGenBuffers)(GLsizei n, GLuint * buffers);
typedef void (*fnptr_glBindBuffer)(GLenum target, GLuint buffer);
typedef void (*fnptr_glBufferData)(GLenum target, GLsizeiptr size, const void * data, GLenum usage);
typedef void (*fnptr_glVertexAttribPointer)(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void * pointer);
typedef void (*fnptr_glEnableVertexAttribArray)(GLuint index);
typedef GLuint (*fnptr_glCreateShader)(GLenum shaderType);
typedef void (*fnptr_glShaderSource)(GLuint shader, GLsizei count, const GLchar **string, const GLint *length);
typedef void (*fnptr_glCompileShader)(GLuint shader);
typedef void (*fnptr_glGetShaderiv)(GLuint shader, GLenum pname, GLint *params);
typedef void (*fnptr_glGetShaderInfoLog)(GLuint shader, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
typedef GLuint (*fnptr_glCreateProgram)(void);
typedef void (*fnptr_glAttachShader)(GLuint program, GLuint shader);
typedef void (*fnptr_glLinkProgram)(GLuint program);
typedef void (*fnptr_glGetProgramiv)(GLuint program, GLenum pname, GLint *params);
typedef void (*fnptr_glGetProgramInfoLog)(GLuint program, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
typedef void (*fnptr_glDeleteShader)(GLuint shader);
typedef void (*fnptr_glUseProgram)(GLuint program);
typedef void (*fnptr_glBufferSubData)(GLenum target, GLintptr offset, GLsizeiptr size, const void * data);
typedef void (*fnptr_glDrawArrays)(GLenum mode, GLint first, GLsizei count);
typedef void (*fnptr_glDrawArraysInstanced)(GLenum mode, GLint first, GLsizei count, GLsizei instancecount);
typedef void (*fnptr_glDrawElementsInstanced)(GLenum mode, GLsizei count, GLenum type, const void * indices, GLsizei instancecount);
typedef GLint (*fnptr_glGetUniformLocation)(GLuint program, const GLchar *name);
typedef void (*fnptr_glUniformMatrix4fv)(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (*fnptr_glUniformMatrix4fv)(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (*fnptr_glVertexAttribDivisor)(GLuint index, GLuint divisor);
typedef void (*fnptr_glGenTextures)(GLsizei n, GLuint * textures);
typedef void (*fnptr_glBindTexture)(GLenum target, GLuint texture);
typedef void (*fnptr_glTexParameteri)(GLenum target, GLenum pname, GLint param);
typedef void (*fnptr_glTexImage2D)(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const void * data);
typedef void (*fnptr_glTexSubImage2D)(GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum format, GLenum type, const void * pixels);


#define VERTEX_ARRAY_BYTES MB(1)
#define INDEX_ARRAY_BYTES MB(1)

struct OpenGLState
{
    HWND hwnd;
    u32_m screen_width;
    u32_m screen_height;

    HGLRC gl_context;
    GLuint last_gl_error;

    GLuint vertex_array_object;
    GLuint index_buffer_object;
    GLuint vertex_buffer_object_vx;
    GLuint vertex_buffer_object_vy;
    GLuint vertex_buffer_object_vz;
    GLuint vertex_buffer_object_nx;
    GLuint vertex_buffer_object_ny;
    GLuint vertex_buffer_object_nz;
    GLuint instanced_vertex_buffer_object_offset_x;
    GLuint instanced_vertex_buffer_object_offset_y;
    GLuint instanced_vertex_buffer_object_offset_z;
    GLuint shader_program;

    GLuint debug_line_vertex_array_object;
    GLuint debug_line_vertex_buffer_object_vertices;
    GLuint debug_line_instanced_vertex_buffer_object_start_x;
    GLuint debug_line_instanced_vertex_buffer_object_start_y;
    GLuint debug_line_instanced_vertex_buffer_object_start_z;
    GLuint debug_line_instanced_vertex_buffer_object_end_x;
    GLuint debug_line_instanced_vertex_buffer_object_end_y;
    GLuint debug_line_instanced_vertex_buffer_object_end_z;
    GLuint debug_line_instanced_vertex_buffer_object_color_r;
    GLuint debug_line_instanced_vertex_buffer_object_color_g;
    GLuint debug_line_instanced_vertex_buffer_object_color_b;
    GLuint debug_line_shader_program;

    GLuint fullscreen_quad_vertex_array_object;
    GLuint fullscreen_quad_vertex_buffer_object_pos;
    GLuint fullscreen_quad_texture;
    GLuint textured_quad_shader_program;

    fnptr_wglCreateContextAttribsARB wglCreateContextAttribsARB;
    fnptr_glGetError glGetError;
    fnptr_glGenVertexArrays glGenVertexArrays;
    fnptr_glBindVertexArray glBindVertexArray;
    fnptr_glGenBuffers glGenBuffers;
    fnptr_glBindBuffer glBindBuffer;
    fnptr_glBufferData glBufferData;
    fnptr_glVertexAttribPointer glVertexAttribPointer;
    fnptr_glEnableVertexAttribArray glEnableVertexAttribArray;
    fnptr_glCreateShader glCreateShader;
    fnptr_glShaderSource glShaderSource;
    fnptr_glCompileShader glCompileShader;
    fnptr_glGetShaderiv glGetShaderiv;
    fnptr_glGetShaderInfoLog glGetShaderInfoLog;
    fnptr_glCreateProgram glCreateProgram;
    fnptr_glAttachShader glAttachShader;
    fnptr_glLinkProgram glLinkProgram;
    fnptr_glGetProgramiv glGetProgramiv;
    fnptr_glGetProgramInfoLog glGetProgramInfoLog;
    fnptr_glDeleteShader glDeleteShader;
    fnptr_glUseProgram glUseProgram;
    fnptr_glBufferSubData glBufferSubData;
    fnptr_glDrawArrays glDrawArrays;
    fnptr_glDrawArraysInstanced glDrawArraysInstanced;
    fnptr_glDrawElementsInstanced glDrawElementsInstanced;
    fnptr_glGetUniformLocation glGetUniformLocation;
    fnptr_glUniformMatrix4fv glUniformMatrix4fv;
    fnptr_glVertexAttribDivisor glVertexAttribDivisor;
    fnptr_glGenTextures glGenTextures;
    fnptr_glBindTexture glBindTexture;
    fnptr_glTexParameteri glTexParameteri;
    fnptr_glTexImage2D glTexImage2D;
    fnptr_glTexSubImage2D glTexSubImage2D;
};

#define CALL_GL(fn, ...) \
do { \
    g_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = g_opengl_state->glGetError(); \
    g_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
} while(0)

#define CALL_GL_RET(out, type, fn, ...) \
do { \
    type MACRO_r = g_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = g_opengl_state->glGetError(); \
    g_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
    *(out) = MACRO_r; \
} while(0)

MAYBE_UNUSED static f32 s_cube_mesh_vx[] = { -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f };
MAYBE_UNUSED static f32 s_cube_mesh_vy[] = { -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f };
MAYBE_UNUSED static f32 s_cube_mesh_vz[] = { -0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f };
MAYBE_UNUSED static u32 s_cube_mesh_indices[] = 
{
    0, 4, 2, // -X
    2, 4, 6, // -X

    3, 5, 1, // +X
    3, 7, 5, // +X

    0, 1, 4, // -Y
    4, 1, 5, // -Y

    2, 6, 7, // +Y
    3, 2, 7, // +Y

    1, 0, 2, // -Z
    2, 3, 1, // -Z

    4, 5, 6, // +Z
    6, 5, 7, // +Z
};

static GLuint make_shader_program(const char* vertex_source, const char* fragment_source)
{
    GLsizei debug_info_len;
    char debug_info_buf[512];
    GLint shader_compile_success;

    GLuint vertex_shader;
    CALL_GL_RET(&vertex_shader, GLuint, glCreateShader, GL_VERTEX_SHADER);
    CALL_GL(glShaderSource, vertex_shader, 1, &vertex_source, NULL);
    CALL_GL(glCompileShader, vertex_shader);
    CALL_GL(glGetShaderiv, vertex_shader, GL_COMPILE_STATUS, &shader_compile_success);
    if(!shader_compile_success)
    {
        CALL_GL(glGetShaderInfoLog, vertex_shader, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
        WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to compile vertex shader.");
    }

    GLuint fragment_shader;
    CALL_GL_RET(&fragment_shader, GLuint, glCreateShader, GL_FRAGMENT_SHADER);
    CALL_GL(glShaderSource, fragment_shader, 1, &fragment_source, NULL);
    CALL_GL(glCompileShader, fragment_shader);
    CALL_GL(glGetShaderiv, fragment_shader, GL_COMPILE_STATUS, &shader_compile_success);
    if(!shader_compile_success)
    {
        CALL_GL(glGetShaderInfoLog, vertex_shader, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
        WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to compile fragment shader.");
    }

    u32_m result_program;
    CALL_GL_RET(&result_program, GLuint, glCreateProgram);
    CALL_GL(glAttachShader, result_program, vertex_shader);
    CALL_GL(glAttachShader, result_program, fragment_shader);
    CALL_GL(glLinkProgram, result_program);
    
    CALL_GL(glGetProgramiv, result_program, GL_LINK_STATUS, &shader_compile_success);
    if(!shader_compile_success)
    {
        CALL_GL(glGetProgramInfoLog, result_program, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
        WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to link shader.");
    }
    
    CALL_GL(glDeleteShader, vertex_shader);
    CALL_GL(glDeleteShader, fragment_shader); 

    return result_program;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Draw array

#define MAX_ARROWS 1024
struct DrawData
{
    u32_m frame_buffer_width;
    u32_m frame_buffer_height;
    u32_m* frame_buffer;

    u32_m num_arrows;
    f32_m arrows_start_x[MAX_ARROWS];
    f32_m arrows_start_y[MAX_ARROWS];
    f32_m arrows_start_z[MAX_ARROWS];
    f32_m arrows_end_x[MAX_ARROWS];
    f32_m arrows_end_y[MAX_ARROWS];
    f32_m arrows_end_z[MAX_ARROWS];
    f32_m arrows_color_r[MAX_ARROWS];
    f32_m arrows_color_g[MAX_ARROWS];
    f32_m arrows_color_b[MAX_ARROWS];
};

MAYBE_UNUSED static inline void init_draw_data(struct DrawData* draw_data, u32 screen_width, u32 screen_height, struct MemoryArena* memory_arena)
{
    draw_data->num_arrows = 0;

    draw_data->frame_buffer_width = screen_width;
    draw_data->frame_buffer_height = screen_height;
    draw_data->frame_buffer = (u32_m*)MEMORY_ARENA_ALLOCATE_ZEROED(memory_arena, screen_width * screen_height * sizeof(u32));
}

MAYBE_UNUSED static inline void reset_draw_data(struct DrawData* draw_data)
{
    draw_data->num_arrows = 0;
    memset(draw_data->frame_buffer, 0, draw_data->frame_buffer_width * draw_data->frame_buffer_height * sizeof(u32));
}

MAYBE_UNUSED static inline void draw_arrow(struct DrawData* draw_data, v3 start, v3 end, v3 color)
{
    u64 n = draw_data->num_arrows++;
    draw_data->arrows_start_x[n] = start.x;
    draw_data->arrows_start_y[n] = start.y;
    draw_data->arrows_start_z[n] = start.z;
    draw_data->arrows_end_x[n] = end.x;
    draw_data->arrows_end_y[n] = end.y;
    draw_data->arrows_end_z[n] = end.z;
    draw_data->arrows_color_r[n] = color.x;
    draw_data->arrows_color_g[n] = color.y;
    draw_data->arrows_color_b[n] = color.z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Input

#define MAX_KEYBOARD_KEYS 256
#define MAX_MOUSE_KEYS 32
struct InputState
{
    u32_m key[MAX_KEYBOARD_KEYS];
    u32_m mouse_key[MAX_MOUSE_KEYS];

    u32_m last_key[MAX_KEYBOARD_KEYS];
    u32_m last_mouse_key[MAX_MOUSE_KEYS];

    u32_m fps_mode; 
    s32_m mouse_screen_pos_x;
    s32_m mouse_screen_pos_y;
    s32_m mouse_screen_dx;
    s32_m mouse_screen_dy;
};

enum KeyboardKey
{
    // https://learn.microsoft.com/en-us/windows/win32/inputdev/virtual-key-codes
    KB_ESCAPE = VK_ESCAPE,
    KB_SPACE = VK_SPACE,
    KB_LCTRL = VK_CONTROL,
    KB_0 = 0x30,    
    KB_1 = 0x31,    
    KB_2 = 0x32,    
    KB_3 = 0x33,    
    KB_4 = 0x34,    
    KB_5 = 0x35,    
    KB_6 = 0x36,    
    KB_7 = 0x37,    
    KB_8 = 0x38,    
    KB_9 = 0x39,    
    KB_A = 0x41,    
    KB_B = 0x42,    
    KB_C = 0x43,    
    KB_D = 0x44,    
    KB_E = 0x45,    
    KB_F = 0x46,    
    KB_G = 0x47,    
    KB_H = 0x48,    
    KB_I = 0x49,    
    KB_J = 0x4A,    
    KB_K = 0x4B,    
    KB_L = 0x4C,    
    KB_M = 0x4D,    
    KB_N = 0x4E,    
    KB_O = 0x4F,    
    KB_P = 0x50,    
    KB_Q = 0x51,    
    KB_R = 0x52,    
    KB_S = 0x53,    
    KB_T = 0x54,    
    KB_U = 0x55,    
    KB_V = 0x56,    
    KB_W = 0x57,    
    KB_X = 0x58,    
    KB_Y = 0x59,    
    KB_Z = 0x5A,    
};

static inline u32 is_key_down(const struct InputState* input_state, enum KeyboardKey k)
{
    ASSERT(k >= 0, "is_key_down key code underflow");
    ASSERT(k < MAX_KEYBOARD_KEYS, "is_key_down key code overflow");
    return input_state->key[(u32)k];
}

static inline u32 is_key_toggled_down(const struct InputState* input_state, enum KeyboardKey k)
{
    ASSERT(k >= 0, "is_key_down key code underflow");
    ASSERT(k < MAX_KEYBOARD_KEYS, "is_key_down key code overflow");
    return input_state->key[(u32)k] && !input_state->last_key[(u32)k];
}

static inline void show_cursor(u32_m show)
{
    if(show)
    {
        // https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-showcursor
        s32_m show_cursor_display_counter = ShowCursor(0);
        while(show_cursor_display_counter >= 0)
        {
            show_cursor_display_counter = ShowCursor(0);
        }
    }
    else
    {
        // https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-showcursor
        s32_m show_cursor_display_counter = ShowCursor(1);
        while(show_cursor_display_counter < 0)
        {
            show_cursor_display_counter = ShowCursor(1);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Time
static s64 get_timestamp_us()
{
    // https://learn.microsoft.com/en-us/windows/win32/sysinfo/acquiring-high-resolution-time-stamps
    LARGE_INTEGER cy;
    QueryPerformanceCounter(&cy);
    cy.QuadPart *= 1000000LL;
    cy.QuadPart /= g_clock_freq;
    s64 result_us = cy.QuadPart;
    ASSERT(result_us >= 0, "get_timestamp_us overflow.");
    return result_us;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



LRESULT CALLBACK WindowProc(HWND window, UINT message, WPARAM wParam, LPARAM lParam)
{
    LRESULT result = 0; 

    switch (message)
    {
        case WM_KEYDOWN: 
        {
            g_input_state->key[wParam] = 1;
            break;
        }

        case WM_KEYUP:
        {
            g_input_state->key[wParam] = 0;
            break;
        }

        case WM_LBUTTONDOWN:
        {
            g_input_state->mouse_key[0] = 1;
            break;
        }

        case WM_LBUTTONUP:
        {
            g_input_state->mouse_key[0] = 0;
            break;
        }

        case WM_RBUTTONDOWN:
        {
            g_input_state->mouse_key[1] = 1;
            break;
        }

        case WM_RBUTTONUP:
        {
            g_input_state->mouse_key[1] = 0;
            break;
        }

        case WM_MOUSEMOVE: 
        {
            struct InputState* input_state = g_input_state;
            if(input_state->fps_mode)
            {
                RECT clip_rect;
                const BOOL get_client_rect_success = GetWindowRect(g_opengl_state->hwnd, &clip_rect);
                ASSERT(get_client_rect_success, "GetWindowRect failed.");
                const BOOL clip_cursor_success = ClipCursor(&clip_rect);
                ASSERT(clip_cursor_success, "ClipCursor failed.");
            }
            else
            {
                const BOOL clip_cursor_success = ClipCursor(NULL);
                ASSERT(clip_cursor_success, "ClipCursor failed.");
            }
            break;
        }

        case WM_SIZE:
        {
            const UINT width = LOWORD(lParam);
            const UINT height = HIWORD(lParam);
            g_opengl_state->screen_width = width;
            g_opengl_state->screen_height = height;
            glViewport(0, 0, width, height);
            break;
        }

        case WM_KILLFOCUS:
        {
            struct InputState* input_state = g_input_state;
            input_state->fps_mode = 0;
            break;
        }

        case WM_DESTROY:
        {
            // https://learn.microsoft.com/en-us/windows/win32/learnwin32/closing-the-window?redirectedfrom=MSDN
            PostQuitMessage(0);
            return 0;
            
        }

        default:
        {
            result = DefWindowProc(window, message, wParam, lParam);
            break;
        }
    }

    return result;
}

static void do_one_frame(struct MemoryArena*);

int WinMainCRTStartup()
{
    // Init logging file.
    g_log_file = CreateFileA("log.txt", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    ASSERT(g_log_file != INVALID_HANDLE_VALUE, "Couldn't open log file for writing.");

    
    // Init application memory.
    void* main_memory_arena_starting_address = (void*)0x100000;
    u64 main_memory_arena_cap = MB(15);
    void* main_memory_arena_storage = VirtualAlloc(
        main_memory_arena_starting_address,
        main_memory_arena_cap,
        MEM_RESERVE | MEM_COMMIT,
        PAGE_READWRITE
    );
    ASSERT(main_memory_arena_storage, "Could not allocate main memory arena.");
    struct MemoryArena main_memory_arena = memory_arena_init(main_memory_arena_storage, main_memory_arena_cap);


    // Init engine graphics state.
    g_opengl_state = (struct OpenGLState*)MEMORY_ARENA_ALLOCATE_ZEROED(&main_memory_arena, sizeof(*g_opengl_state));
    {
        const HMODULE hinstance = GetModuleHandle(0);

        WNDCLASS window_class = {0};
        window_class.style = CS_HREDRAW|CS_VREDRAW|CS_OWNDC;
        window_class.lpfnWndProc = WindowProc;
        window_class.hInstance = hinstance;
        window_class.hCursor = LoadCursor(NULL, IDC_ARROW);
        window_class.lpszClassName = "Windows Program Class";
        ATOM register_class_result = RegisterClass(&window_class);
        ASSERT(register_class_result, "Register class failed");

        u64 window_width = 1920;
        u64 window_height = 1080;
        u64 window_start_x = (GetSystemMetrics(SM_CXSCREEN) - window_width) / 2;
        u64 window_start_y = (GetSystemMetrics(SM_CYSCREEN) - window_height) / 2;
        g_opengl_state->hwnd = CreateWindowEx(0,          // Extended style
                window_class.lpszClassName,               // Class name
                "",                                       // Window name
                //WS_OVERLAPPEDWINDOW | WS_VISIBLE,       // Style of the window
                WS_POPUP | WS_VISIBLE,                    // Style of the window
                (int)window_start_x,                      // Initial X position
                (int)window_start_y,                      // Initial Y position
                (int)window_width,                        // Initial width
                (int)window_height,                       // Initial height
                0,                                        // Handle to the window parent
                0,                                        // Handle to a menu
                hinstance,                                // Handle to an instance
                0);
        ASSERT(g_opengl_state->hwnd != NULL, "Failed to create a window");

        const HDC dc = GetDC(g_opengl_state->hwnd);
        ASSERT(dc != NULL, "GetDC failed.");

        RECT client_rect;
        const BOOL get_client_rect_success = GetClientRect(g_opengl_state->hwnd, &client_rect);
        ASSERT(get_client_rect_success, "GetClientRect failed.");
        g_opengl_state->screen_width = client_rect.right;
        g_opengl_state->screen_height = client_rect.bottom;

        BOOL clip_cursor_success = ClipCursor(&client_rect);
        ASSERT(clip_cursor_success, "ClipCursor failed.");

        HMODULE opengl32_dll_module = LoadLibraryA("opengl32.dll");
        ASSERT(opengl32_dll_module != NULL, "Could not load opengl32.dll");

        PIXELFORMATDESCRIPTOR pfd = {0};
        pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
        pfd.nVersion = 1;
        pfd.dwFlags = PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW;
        pfd.iPixelType = PFD_TYPE_RGBA;
        pfd.cColorBits = 32;
        pfd.cDepthBits = 32;
        pfd.iLayerType = PFD_MAIN_PLANE;
        const int pixel_format = ChoosePixelFormat(dc, &pfd);
        ASSERT(pixel_format, "ChoosePixelFormat failure");

        const BOOL set_pixel_format_success = SetPixelFormat(dc, pixel_format, &pfd);
        ASSERT(set_pixel_format_success, "SetPixelFormat failure");

        HGLRC temp_context = wglCreateContext(dc);
        ASSERT(temp_context != NULL, "Could not make opengl context.");

        const BOOL wgl_make_current_temp_success = wglMakeCurrent(dc, temp_context);
        ASSERT(wgl_make_current_temp_success == TRUE, "wglMakeCurrent failed.");

        g_opengl_state->wglCreateContextAttribsARB = (fnptr_wglCreateContextAttribsARB)load_gl_fn(opengl32_dll_module, "wglCreateContextAttribsARB");

        int attribs[] =
        {
            WGL_CONTEXT_MAJOR_VERSION_ARB, 3,
            WGL_CONTEXT_MINOR_VERSION_ARB, 3,
            WGL_CONTEXT_FLAGS_ARB, 0,
            0
        };
        g_opengl_state->gl_context = g_opengl_state->wglCreateContextAttribsARB(dc, 0, attribs);
        ASSERT(g_opengl_state->gl_context != NULL, "Could not make opengl context.");

        const BOOL wgl_make_current_null_success = wglMakeCurrent(NULL, NULL);
        ASSERT(wgl_make_current_null_success == TRUE, "wglMakeCurrent failed.");

        const BOOL wgl_make_currentxt_success = wglDeleteContext(temp_context);
        ASSERT(wgl_make_currentxt_success == TRUE, "wglDeleteContext failed.");

        const BOOL wgl_make_current_final_success = wglMakeCurrent(dc, g_opengl_state->gl_context);
        ASSERT(wgl_make_current_final_success == TRUE, "wglMakeCurrent failed.");

        // Load OpenGL functions.
        g_opengl_state->glGetError = (fnptr_glGetError)load_gl_fn(opengl32_dll_module, "glGetError");
        g_opengl_state->glGenVertexArrays = (fnptr_glGenVertexArrays)load_gl_fn(opengl32_dll_module, "glGenVertexArrays");
        g_opengl_state->glBindVertexArray = (fnptr_glBindVertexArray)load_gl_fn(opengl32_dll_module, "glBindVertexArray");
        g_opengl_state->glGenBuffers = (fnptr_glGenBuffers)load_gl_fn(opengl32_dll_module, "glGenBuffers");
        g_opengl_state->glBindBuffer = (fnptr_glBindBuffer)load_gl_fn(opengl32_dll_module, "glBindBuffer");
        g_opengl_state->glBufferData = (fnptr_glBufferData)load_gl_fn(opengl32_dll_module, "glBufferData");
        g_opengl_state->glVertexAttribPointer = (fnptr_glVertexAttribPointer)load_gl_fn(opengl32_dll_module, "glVertexAttribPointer");
        g_opengl_state->glEnableVertexAttribArray = (fnptr_glEnableVertexAttribArray)load_gl_fn(opengl32_dll_module, "glEnableVertexAttribArray");
        g_opengl_state->glCreateShader = (fnptr_glCreateShader)load_gl_fn(opengl32_dll_module, "glCreateShader");
        g_opengl_state->glShaderSource = (fnptr_glShaderSource)load_gl_fn(opengl32_dll_module, "glShaderSource");
        g_opengl_state->glCompileShader = (fnptr_glCompileShader)load_gl_fn(opengl32_dll_module, "glCompileShader");
        g_opengl_state->glGetShaderiv = (fnptr_glGetShaderiv)load_gl_fn(opengl32_dll_module, "glGetShaderiv");
        g_opengl_state->glGetShaderInfoLog = (fnptr_glGetShaderInfoLog)load_gl_fn(opengl32_dll_module, "glGetShaderInfoLog");
        g_opengl_state->glCreateProgram = (fnptr_glCreateProgram)load_gl_fn(opengl32_dll_module, "glCreateProgram");
        g_opengl_state->glAttachShader = (fnptr_glAttachShader)load_gl_fn(opengl32_dll_module, "glAttachShader");
        g_opengl_state->glLinkProgram = (fnptr_glLinkProgram)load_gl_fn(opengl32_dll_module, "glLinkProgram");
        g_opengl_state->glGetProgramiv = (fnptr_glGetProgramiv)load_gl_fn(opengl32_dll_module, "glGetProgramiv");
        g_opengl_state->glGetProgramInfoLog = (fnptr_glGetProgramInfoLog)load_gl_fn(opengl32_dll_module, "glGetProgramInfoLog");
        g_opengl_state->glDeleteShader = (fnptr_glDeleteShader)load_gl_fn(opengl32_dll_module, "glDeleteShader");
        g_opengl_state->glUseProgram = (fnptr_glUseProgram)load_gl_fn(opengl32_dll_module, "glUseProgram");
        g_opengl_state->glBufferSubData = (fnptr_glBufferSubData)load_gl_fn(opengl32_dll_module, "glBufferSubData");
        g_opengl_state->glDrawArrays = (fnptr_glDrawArrays)load_gl_fn(opengl32_dll_module, "glDrawArrays");
        g_opengl_state->glDrawArraysInstanced = (fnptr_glDrawArraysInstanced)load_gl_fn(opengl32_dll_module, "glDrawArraysInstanced");
        g_opengl_state->glDrawElementsInstanced = (fnptr_glDrawElementsInstanced)load_gl_fn(opengl32_dll_module, "glDrawElementsInstanced");
        g_opengl_state->glGetUniformLocation = (fnptr_glGetUniformLocation)load_gl_fn(opengl32_dll_module, "glGetUniformLocation");
        g_opengl_state->glUniformMatrix4fv = (fnptr_glUniformMatrix4fv)load_gl_fn(opengl32_dll_module, "glUniformMatrix4fv");
        g_opengl_state->glVertexAttribDivisor = (fnptr_glVertexAttribDivisor)load_gl_fn(opengl32_dll_module, "glVertexAttribDivisor");

        g_opengl_state->glGenTextures = (fnptr_glGenTextures)load_gl_fn(opengl32_dll_module, "glGenTextures");
        g_opengl_state->glBindTexture = (fnptr_glBindTexture)load_gl_fn(opengl32_dll_module, "glBindTexture");
        g_opengl_state->glTexParameteri = (fnptr_glTexParameteri)load_gl_fn(opengl32_dll_module, "glTexParameteri");
        g_opengl_state->glTexImage2D = (fnptr_glTexImage2D)load_gl_fn(opengl32_dll_module, "glTexImage2D");
        g_opengl_state->glTexSubImage2D = (fnptr_glTexSubImage2D)load_gl_fn(opengl32_dll_module, "glTexSubImage2D");

        glViewport(0, 0, (GLsizei)g_opengl_state->screen_width, (GLsizei)g_opengl_state->screen_height);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CCW);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);


        {
            // Vertex array object
            CALL_GL(glGenVertexArrays, 1, &g_opengl_state->vertex_array_object);
            CALL_GL(glBindVertexArray, g_opengl_state->vertex_array_object);

            // Index buffer object
            CALL_GL(glGenBuffers, 1, &g_opengl_state->index_buffer_object);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, g_opengl_state->index_buffer_object);
            CALL_GL(glBufferData, GL_ELEMENT_ARRAY_BUFFER, INDEX_ARRAY_BYTES, 0, GL_DYNAMIC_DRAW);

            u32_m attr_idx = 0;

            // Vertex buffer: vx
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_vx);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vx);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: vy
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_vy);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vy);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: vz
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_vz);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vz);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: nx
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_nx);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_nx);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: ny
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_ny);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_ny);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: nz
            CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_nz);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_nz);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Instanced vertex buffer: offset X
            CALL_GL(glGenBuffers, 1, &g_opengl_state->instanced_vertex_buffer_object_offset_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Instanced vertex buffer: offset Y
            CALL_GL(glGenBuffers, 1, &g_opengl_state->instanced_vertex_buffer_object_offset_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Instanced vertex buffer: offset Z
            CALL_GL(glGenBuffers, 1, &g_opengl_state->instanced_vertex_buffer_object_offset_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

            // Shaders
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in float a_vx;\n"
                "layout (location = 1) in float a_vy;\n"
                "layout (location = 2) in float a_vz;\n"
                "layout (location = 3) in float a_nx;\n"
                "layout (location = 4) in float a_ny;\n"
                "layout (location = 5) in float a_nz;\n"
                "layout (location = 6) in float a_offset_x;\n"
                "layout (location = 7) in float a_offset_y;\n"
                "layout (location = 8) in float a_offset_z;\n"
                "uniform mat4 m_mvp;\n"
                "out vec3 v_color;\n"
                "void main()\n"
                "{\n"
                "    vec3 pos = vec3(a_vx, a_vy, a_vz) + vec3(a_offset_x, a_offset_y, a_offset_z);\n"
                "    v_color = vec3(a_nx, a_ny, a_nz);"
                "    gl_Position = m_mvp * vec4(pos, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "in vec3 v_color;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    result_frag_color = vec4(v_color, 1.0f);\n"
                //"    result_frag_color = vec4(0.05f, 0.5f, 0.0f, 1.0f);\n"
                //"    result_frag_color = vec4(0.05f, 0.5f, 0.0f, 1.0f);\n"
                "}\n"
                "\n";
            g_opengl_state->shader_program = make_shader_program(vertex_source, fragment_source);
        }


        {
            // Debug lines
            // Vertex array object
            CALL_GL(glGenVertexArrays, 1, &g_opengl_state->debug_line_vertex_array_object);
            CALL_GL(glBindVertexArray, g_opengl_state->debug_line_vertex_array_object);

            u32_m attr_idx = 0;

            // Vertex buffer: vertices
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_vertex_buffer_object_vertices);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_vertex_buffer_object_vertices);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, 6 * sizeof(f32), NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 3, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: start x
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: start y
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: start z
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end x
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end y
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end z
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color r
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color g
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color b
            CALL_GL(glGenBuffers, 1, &g_opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

            // Shaders
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in vec3 a_vertex;\n"
                "layout (location = 1) in float a_start_x;\n"
                "layout (location = 2) in float a_start_y;\n"
                "layout (location = 3) in float a_start_z;\n"
                "layout (location = 4) in float a_end_x;\n"
                "layout (location = 5) in float a_end_y;\n"
                "layout (location = 6) in float a_end_z;\n"
                "layout (location = 7) in float a_color_r;\n"
                "layout (location = 8) in float a_color_g;\n"
                "layout (location = 9) in float a_color_b;\n"
                "uniform mat4 m_mvp;\n"
                "out vec3 v_color;\n"
                "void main()\n"
                "{\n"
                "    vec3 start = vec3(a_start_x, a_start_y, a_start_z);\n"
                "    vec3 end = vec3(a_end_x, a_end_y, a_end_z);\n"
                "    vec3 pos = mix(start, end, a_vertex);\n"
                "    v_color = vec3(a_color_r, a_color_g, a_color_b);\n"
                "    gl_Position = m_mvp * vec4(pos, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "in vec3 v_color;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    result_frag_color = vec4(v_color, 1.0f);\n"
                "}\n"
                "\n";
            g_opengl_state->debug_line_shader_program = make_shader_program(vertex_source, fragment_source);
        }


        {
            CALL_GL(glGenTextures, 1, &g_opengl_state->fullscreen_quad_texture);
            CALL_GL(glBindTexture, GL_TEXTURE_2D, g_opengl_state->fullscreen_quad_texture);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            CALL_GL(glTexImage2D, GL_TEXTURE_2D,
                         0,
                         GL_RGBA,
                         g_opengl_state->screen_width,
                         g_opengl_state->screen_height,
                         0,
                         GL_RGBA,
                         GL_UNSIGNED_BYTE,
                         NULL);
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in vec2 a_pos;\n"
                "out vec2 v_uv;\n"
                "void main()\n"
                "{\n"
                "    v_uv = a_pos * 0.5f + vec2(0.5f, 0.5f);\n"
                "    gl_Position = vec4(a_pos, 0.0f, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "uniform sampler2D tex_sampler;\n"
                "in vec2 v_uv;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    result_frag_color = texture(tex_sampler, v_uv);\n"
                "}\n"
                "\n";
            g_opengl_state->textured_quad_shader_program = make_shader_program(vertex_source, fragment_source);

            CALL_GL(glGenVertexArrays, 1, &g_opengl_state->fullscreen_quad_vertex_array_object);
            CALL_GL(glBindVertexArray, g_opengl_state->fullscreen_quad_vertex_array_object);

            u32_m attr_idx = 0;

            // Vertex buffer: vertices
            CALL_GL(glGenBuffers, 1, &g_opengl_state->fullscreen_quad_vertex_buffer_object_pos);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->fullscreen_quad_vertex_buffer_object_pos);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, 6 * 2 * sizeof(f32), NULL, GL_STATIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 2, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
        }
    }


    // Init draw array.
    g_draw_data = (struct DrawData*)MEMORY_ARENA_ALLOCATE(&main_memory_arena, sizeof(struct DrawData));
    init_draw_data(g_draw_data, g_opengl_state->screen_width, g_opengl_state->screen_height, &main_memory_arena);
    reset_draw_data(g_draw_data);


    // Init engine input state.
    // Dependent on window being created.
    {
        g_input_state = (struct InputState*)MEMORY_ARENA_ALLOCATE_ZEROED(&main_memory_arena, sizeof(*g_input_state));
        struct InputState* input_state = g_input_state;

        POINT p;
        const BOOL get_cursor_pos_succes = GetCursorPos(&p);
        ASSERT(get_cursor_pos_succes, "Failed GetCursorPos");
        const BOOL screen_to_client_success = ScreenToClient(g_opengl_state->hwnd, &p);
        ASSERT(screen_to_client_success, "Failed ScreenToClient");
        input_state->mouse_screen_pos_x = (s32_m)p.x;
        input_state->mouse_screen_pos_y = (s32_m)p.y;
        input_state->fps_mode = 0;
    }

    // Init engine time.
    {
        LARGE_INTEGER clock_freq;
        QueryPerformanceFrequency(&clock_freq);
        g_clock_freq = clock_freq.QuadPart;
    }


    u64 frame_memory_arena_bytes = MB(5);
    struct MemoryArena frame_memory_arena = memory_arena_init(MEMORY_ARENA_ALLOCATE(&main_memory_arena, frame_memory_arena_bytes), frame_memory_arena_bytes);


    // Main loop
    s64 engine_start_timestamp_us = get_timestamp_us();
    s64_m last_frame_timestamp_us = engine_start_timestamp_us;
    while(1)
    {


        MSG msg;
        while(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
        {
            // https://learn.microsoft.com/en-us/windows/win32/learnwin32/closing-the-window?redirectedfrom=MSDN
            if(msg.message == WM_QUIT)
            {
                ExitProcess(0);
            }
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        if(is_key_down(g_input_state, KB_ESCAPE))
        {
            ExitProcess(0);
        }


        // Main loop time control.
        {
            s64 time_since_last_frame_us = get_timestamp_us() - last_frame_timestamp_us;
            if(time_since_last_frame_us >= ENGINE_FRAME_DURATION_US)
            {
                last_frame_timestamp_us += ENGINE_FRAME_DURATION_US;
            }
            else
            {
                continue;
            }
        }

        do_one_frame(&frame_memory_arena);

        memory_arena_reset(&frame_memory_arena);

        const HDC dc = GetDC(g_opengl_state->hwnd);
        BOOL swap_buffers_success = SwapBuffers(dc);
        ASSERT(swap_buffers_success, "SwapBuffers failed.");
    }

    CloseHandle(g_log_file);

    ExitProcess(0);
}


static void do_one_frame(struct MemoryArena* memory_arena)
{
    {
        struct InputState* input_state = g_input_state;

        if(is_key_toggled_down(input_state, KB_E))
        {
            input_state->fps_mode = !input_state->fps_mode;

            RECT clip_rect;
            const BOOL get_client_rect_success = GetClientRect(g_opengl_state->hwnd, &clip_rect);
            ASSERT(get_client_rect_success, "GetClientRect failed.");
            s32 tx = (clip_rect.left + clip_rect.right) / 2;
            s32 ty = (clip_rect.bottom + clip_rect.top) / 2;
            POINT client_to_screen_point = {
                .x = tx,
                .y = ty,
            };
            const BOOL set_cursor_pos_success = SetCursorPos(client_to_screen_point.x, client_to_screen_point.y);
            ASSERT(set_cursor_pos_success, "SetCursorPos failed.");
            input_state->mouse_screen_pos_x = client_to_screen_point.x;
            input_state->mouse_screen_pos_y = client_to_screen_point.y;
        }


        s32 last_mouse_screen_pos_x = input_state->mouse_screen_pos_x;
        s32 last_mouse_screen_pos_y = input_state->mouse_screen_pos_y;
        POINT p;
        const BOOL get_cursor_pos_succes = GetCursorPos(&p);
        ASSERT(get_cursor_pos_succes, "Failed GetCursorPos");
        input_state->mouse_screen_pos_x = (s32_m)p.x;
        input_state->mouse_screen_pos_y = (s32_m)p.y;

        if(input_state->fps_mode)
        {
            show_cursor(1);
            RECT clip_rect;
            const BOOL get_client_rect_success = GetClientRect(g_opengl_state->hwnd, &clip_rect);
            ASSERT(get_client_rect_success, "GetClientRect failed.");
            s32 tx = (clip_rect.left + clip_rect.right) / 2;
            s32 ty = (clip_rect.bottom + clip_rect.top) / 2;

            POINT client_to_screen_point = {
                .x = tx,
                .y = ty,
            };
            const BOOL client_to_screen_success = ClientToScreen(g_opengl_state->hwnd, &client_to_screen_point);
            ASSERT(client_to_screen_success, "ClientToScreen failed.");

            input_state->mouse_screen_dx = input_state->mouse_screen_pos_x - last_mouse_screen_pos_x;
            input_state->mouse_screen_dy = input_state->mouse_screen_pos_y - last_mouse_screen_pos_y;

            const BOOL set_cursor_pos_success = SetCursorPos(client_to_screen_point.x, client_to_screen_point.y);
            ASSERT(set_cursor_pos_success, "SetCursorPos failed.");

            input_state->mouse_screen_pos_x = client_to_screen_point.x;
            input_state->mouse_screen_pos_y = client_to_screen_point.y;
        }
        else
        {
            show_cursor(0);
            input_state->mouse_screen_dx = input_state->mouse_screen_pos_x - last_mouse_screen_pos_x;
            input_state->mouse_screen_dy = input_state->mouse_screen_pos_y - last_mouse_screen_pos_y;
        }
    }


    reset_draw_data(g_draw_data);


    // Draw world basis.
    {
        v3 s = {.m={0.0f, 0.0f, 0.0f}};
        v3_m e = {.m={1.0f, 0.0f, 0.0f}};
        draw_arrow(g_draw_data, s, e, e);
        e.x = 0.0f;
        e.y = 1.0f;
        e.z = 0.0f;
        draw_arrow(g_draw_data, s, e, e);
        e.x = 0.0f;
        e.y = 0.0f;
        e.z = 1.0f;
        draw_arrow(g_draw_data, s, e, e);
    }



    u32 terrain_width = 128;
    u32 terrain_vertex_stride = terrain_width + 1;
    u32 num_terrain_vertices = (terrain_width + 1) * (terrain_width + 1);
    f32_m* terrain_vx = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    f32_m* terrain_vy = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    f32_m* terrain_vz = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    f32_m* terrain_nx = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    f32_m* terrain_ny = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    f32_m* terrain_nz = (f32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, num_terrain_vertices * sizeof(f32));
    u32_m* terrain_indices = (u32_m*)MEMORY_ARENA_ALLOCATE(memory_arena, terrain_width * terrain_width * 6 * sizeof(u32));

    u32_m num_terrain_indices = 0;
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
                    terrain_vx[i_z * terrain_vertex_stride + i_x + lane] = x + (f32)lane;
                    terrain_vy[i_z * terrain_vertex_stride + i_x + lane] = n * 100.0f;
                    terrain_vz[i_z * terrain_vertex_stride + i_x + lane] = z;

                    u32 c = rand_u32((u32)(i_z * terrain_vertex_stride + i_x + lane));
                    terrain_nx[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >> 16) & 0xFF) / 256.0f;
                    terrain_ny[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >>  8) & 0xFF) / 256.0f;
                    terrain_nz[i_z * terrain_vertex_stride + i_x + lane] = (f32)((c >>  0) & 0xFF) / 256.0f;
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
                    terrain_vx[i_z * terrain_vertex_stride + i_x] = x + (f32)lane;
                    terrain_vy[i_z * terrain_vertex_stride + i_x] = n * 100.0f;
                    terrain_vz[i_z * terrain_vertex_stride + i_x] = z;

                    terrain_nx[i_z * terrain_vertex_stride + i_x] = n;
                    terrain_ny[i_z * terrain_vertex_stride + i_x] = n;
                    terrain_nz[i_z * terrain_vertex_stride + i_x] = n;
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

                terrain_indices[num_terrain_indices++] = bl;
                terrain_indices[num_terrain_indices++] = br;
                terrain_indices[num_terrain_indices++] = tl;

                terrain_indices[num_terrain_indices++] = tl;
                terrain_indices[num_terrain_indices++] = br;
                terrain_indices[num_terrain_indices++] = tr;
            }
        }
    }


    f32 offset_x[] = { 0.0f };
    f32 offset_y[] = { -1.0f };
    f32 offset_z[] = { 0.0f };

    glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



    {
        // Camera controls
        static f32_m cam_pitch_turns = 0.0f;
        cam_pitch_turns += (float)is_key_down(g_input_state, KB_I) * 0.001f;
        cam_pitch_turns -= (float)is_key_down(g_input_state, KB_K) * 0.001f;


        static f32_m cam_yaw_turns = 0.0f;
        cam_yaw_turns += (float)is_key_down(g_input_state, KB_J) * 0.001f;
        cam_yaw_turns -= (float)is_key_down(g_input_state, KB_L) * 0.001f;

        if(g_input_state->fps_mode)
        {
            cam_pitch_turns -= (f32)(g_input_state->mouse_screen_dy) * 0.0005f;
            cam_yaw_turns -= (f32)(g_input_state->mouse_screen_dx) * 0.0005f;
        }

        static v3_m cam_pos = {.m = {0.0f, 1.0f, 3.0f}};

        mtx4x4_m y_rot_mtx;
        make_y_axis_rotation_mtx(&y_rot_mtx, -cam_yaw_turns);

        mtx4x4_m x_rot_mtx;
        make_x_axis_rotation_mtx(&x_rot_mtx, -cam_pitch_turns);

        mtx4x4_m translation_mtx;
        make_translation_mtx(&translation_mtx, v3_scale(cam_pos, -1.0f));

        mtx4x4_m world_to_cam_mtx_temp;
        mtx4x4_mul(&world_to_cam_mtx_temp, &y_rot_mtx, &translation_mtx);

        mtx4x4_m world_to_cam_mtx;
        mtx4x4_mul(&world_to_cam_mtx, &x_rot_mtx, &world_to_cam_mtx_temp);

        v3 I_cam = {.m={world_to_cam_mtx.m[0*4 + 0], world_to_cam_mtx.m[0*4 + 1], world_to_cam_mtx.m[0*4 + 2]}};
        v3 J_cam = {.m={world_to_cam_mtx.m[1*4 + 0], world_to_cam_mtx.m[1*4 + 1], world_to_cam_mtx.m[1*4 + 2]}};
        v3 K_cam = {.m={world_to_cam_mtx.m[2*4 + 0], world_to_cam_mtx.m[2*4 + 1], world_to_cam_mtx.m[2*4 + 2]}};

        f32 speed = 0.5f;
        cam_pos = v3_add(cam_pos, v3_scale(K_cam, -speed * (float)is_key_down(g_input_state, KB_W)));
        cam_pos = v3_add(cam_pos, v3_scale(K_cam,  speed * (float)is_key_down(g_input_state, KB_S)));

        cam_pos = v3_add(cam_pos, v3_scale(I_cam,  speed * (float)is_key_down(g_input_state, KB_D)));
        cam_pos = v3_add(cam_pos, v3_scale(I_cam, -speed * (float)is_key_down(g_input_state, KB_A)));
        
        cam_pos = v3_add(cam_pos, v3_scale(J_cam,  speed * (float)is_key_down(g_input_state, KB_SPACE)));
        cam_pos = v3_add(cam_pos, v3_scale(J_cam, -speed * (float)is_key_down(g_input_state, KB_LCTRL)));

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

        f32 aspect_ratio = (float)g_opengl_state->screen_height / (float)g_opengl_state->screen_width;

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

        mtx4x4_m mvp_mtx;
        mtx4x4_mul(&mvp_mtx, &proj_mtx, &world_to_cam_mtx);

        CALL_GL(glUseProgram, g_opengl_state->shader_program);

        GLint loc;
        CALL_GL_RET(&loc, GLint, glGetUniformLocation, g_opengl_state->shader_program, "m_mvp");
        CALL_GL(glUniformMatrix4fv, loc, 1, 1, &mvp_mtx.m[0]);
        ASSERT(loc != -1, "Failed to bind uniform.");

        CALL_GL(glBindVertexArray, g_opengl_state->vertex_array_object);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vx);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_vx);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vy);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_vy);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_vz);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_vz);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_nx);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_nx);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_ny);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_ny);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_nz);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_terrain_vertices * sizeof(f32), terrain_nz);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_x);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(offset_x), offset_x);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_y);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(offset_y), offset_y);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->instanced_vertex_buffer_object_offset_z);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(offset_z), offset_z);

        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, g_opengl_state->index_buffer_object);
        CALL_GL(glBufferSubData, GL_ELEMENT_ARRAY_BUFFER, 0, num_terrain_indices * sizeof(u32), terrain_indices);

        CALL_GL(glDrawElementsInstanced, GL_TRIANGLES, num_terrain_indices, GL_UNSIGNED_INT, 0, 1/*batch_size*/);
        //CALL_GL(glDrawElementsInstanced, GL_LINES, num_terrain_indices, GL_UNSIGNED_INT, 0, 1/*batch_size*/);

        CALL_GL(glBindVertexArray, 0);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
        CALL_GL(glUseProgram, 0);


        // Render debug lines.
        {
            CALL_GL(glUseProgram, g_opengl_state->debug_line_shader_program);
            
            GLint debug_line_loc;
            CALL_GL_RET(&debug_line_loc, GLint, glGetUniformLocation, g_opengl_state->shader_program, "m_mvp");
            CALL_GL(glUniformMatrix4fv, loc, 1, 1, &mvp_mtx.m[0]);
            ASSERT(debug_line_loc != -1, "Failed to bind uniform.");
            
            CALL_GL(glBindVertexArray, g_opengl_state->debug_line_vertex_array_object);

            f32 line_vertices[6] = {
                0.0f, 0.0f, 0.0f,
                1.0f, 1.0f, 1.0f,
            };
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_vertex_buffer_object_vertices);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(line_vertices), line_vertices);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_start_x), g_draw_data->arrows_start_x);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_start_y), g_draw_data->arrows_start_y);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_start_z), g_draw_data->arrows_start_z);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_end_x), g_draw_data->arrows_end_x);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_end_y), g_draw_data->arrows_end_y);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_end_z), g_draw_data->arrows_end_z);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_color_r), g_draw_data->arrows_color_r);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_color_g), g_draw_data->arrows_color_g);

            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, g_draw_data->num_arrows * sizeof(*g_draw_data->arrows_color_b), g_draw_data->arrows_color_b);
            
            CALL_GL(glDrawArraysInstanced, GL_LINES, 0, 2, g_draw_data->num_arrows);

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
            CALL_GL(glUseProgram, 0);
        }

        // Debug texture
        {

            // random number betwee 0 and 4294967295
            // 1026793478
            // 2550638353
            // 513730960
            // 980227355
            // 1575412276
            // 668711287
            // 3002072601
            // 1157920987

            CALL_GL(glUseProgram, g_opengl_state->textured_quad_shader_program);

#if 0
            {
                u32 random_seed_8[] = {
                    1026793478,
                    2550638353,
                    513730960,
                    980227355,
                    1575412276,
                    668711287,
                    3002072601,
                    1157920987,
                };
                __m256i vrandom = _mm256_loadu_si256((__m256i*)random_seed_8);
                for(u64_m y = 0; y < g_draw_data->frame_buffer_height; y++)
                {
                    for(u64_m x = 0; x < g_draw_data->frame_buffer_width - 8; x += 8)
                    {
                        u32_m storage8[8];
                        _mm256_storeu_si256((__m256i*)storage8, vrandom);

                        for(u64_m lane = 0; lane < 8; lane++)
                        {
                            u32 b = storage8[lane] & 0xFF;
                            g_draw_data->frame_buffer[y * g_draw_data->frame_buffer_width + x + lane] = 0xFF << 24 | b << 16 | b << 8 | (b);
                        }

                        vrandom = rand8_u32(vrandom);
                    }
                }
            }
#elif 0
            {
                f32 offset = 32.0f;
                f32 scale = g_input_state->mouse_screen_pos_x * 0.003f;
                __m256 accum_x = _mm256_fmadd_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(scale), _mm256_set1_ps(offset));
                __m256 accum_y = _mm256_set1_ps(offset);
                for(u64_m y = 0; y < g_draw_data->frame_buffer_height; y++)
                {
                    for(u64_m x = 0; x < g_draw_data->frame_buffer_width - 8; x += 8)
                    {
                        __m256 v = pnoise8(accum_x, accum_y, _mm256_set1_ps(0.0f));

                        f32_m storage8[8];
                        _mm256_storeu_ps(storage8, v);

                        for(u64_m lane = 0; lane < 8; lane++)
                        {
                            u32 b = (u32)((storage8[lane] * 0.5f + 0.5f) * 255.0f);
                            u32 color = 0xFF << 24 | b << 16 | b << 8 | (b);
                            g_draw_data->frame_buffer[y * g_draw_data->frame_buffer_width + x + lane] = color;
                        }

                        accum_x = _mm256_add_ps(accum_x, _mm256_set1_ps(scale * 8.0f));
                    }
                    accum_y = _mm256_add_ps(accum_y, _mm256_set1_ps(scale));
                    accum_x = _mm256_fmadd_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(scale), _mm256_set1_ps(offset));
                }
            }
#endif
    
#if 0
            {
                for(u64_m y = 0; y < g_draw_data->frame_buffer_height; y++)
                {
                    for(u64_m x = 0; x < g_draw_data->frame_buffer_width; x++)
                    {
                        f32 in = (f32)x * 0.04f;
                        __m256 v = approx_sin8(_mm256_set1_ps(in));
                        f32_m storage8[8];
                        _mm256_storeu_ps(storage8, v);
                        u32 b = (u32)((storage8[0] * 0.5f + 0.5f) * 255.0f);
                        g_draw_data->frame_buffer[y * g_draw_data->frame_buffer_width + x] = 0xFF << 24 | b << 16 | b << 8 | (b);
                    }
                }
            }
#endif



#if 0
            CALL_GL(glTexSubImage2D,
                    GL_TEXTURE_2D, // GLenum target
                    0, // GLint level
                    0, // GLint xoffset
                    0, // GLint yoffset
                    g_draw_data->frame_buffer_width,  // GLsizei width
                    g_draw_data->frame_buffer_height, // GLsizei height
                    GL_RGBA,          // GLenum format
                    GL_UNSIGNED_BYTE, // GLenum type
                    g_draw_data->frame_buffer);
            const GLuint asdferr = g_opengl_state->glGetError();
            ASSERT(asdferr == 0, "failed"); \

            CALL_GL(glBindVertexArray, g_opengl_state->fullscreen_quad_vertex_array_object);

            f32 quad_vertices[] = {
                -1.0f, -1.0f,
                 1.0f, -1.0f,
                 1.0f,  1.0f,

                -1.0f, -1.0f,
                 1.0f,  1.0f,
                -1.0f,  1.0f,
            };
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->fullscreen_quad_vertex_buffer_object_pos);
            CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(quad_vertices), quad_vertices);

            CALL_GL(glDrawArrays, GL_TRIANGLES, 0, 6);

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
            CALL_GL(glUseProgram, 0);
#endif
        }
    }

    {
        struct InputState* input_state = g_input_state;
        memcpy(input_state->last_key, input_state->key, sizeof(input_state->last_key));
        memcpy(input_state->last_mouse_key, input_state->mouse_key, sizeof(input_state->last_mouse_key));
    }
}


