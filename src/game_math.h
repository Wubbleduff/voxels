
#pragma once

#include "common.h"

#include <cstdlib>
#include <cmath>

#pragma warning(disable:4201)
#pragma warning(disable:4505)

#define M_PI 3.14159265f

struct v2
{
    union
    {
        struct { float x, y; };
        float v[2];
    };
    
    v2() : x(0.0f), y(0.0f) { }
    v2(const v2 &v) : x(v.x), y(v.y) { }
    v2(float in_x, float in_y) : x(in_x), y(in_y) { }
    v2(float c[2]) : x(c[0]), y(c[1]) { }
};

struct v3
{
    union
    {
        struct { float x; float y; float z; };
        struct { float r; float g; float b; };
        struct { float h; float s; float v; };
        float m[3];
    };
    v3() : x(0.0f), y(0.0f), z(0.0f) { }
    v3(const v3 &v) : x(v.x), y(v.y), z(v.z) { }
    v3(float in_x, float in_y, float in_z) : x(in_x), y(in_y), z(in_z) { }
    v3(v2 v, float a) : x(v.x), y(v.y), z(a)   { }
    v3(float a, v2 v) : x(a),   y(v.x), z(v.y) { }
};

struct v4
{
    union
    {
        struct { float x, y, z, w; };
        struct { float r, g, b, a; };
        struct { float h, s, v, aa; };
        float m[4];
    };
    v4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) { }
    v4(const v4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
    v4(float in_x, float in_y, float in_z, float in_w) : x(in_x), y(in_y), z(in_z), w(in_w) { }
    v4(v2 v, float a, float b) : x(v.x), y(v.y), z(a),   w(b)   { }
    v4(float a, v2 v, float b) : x(a),   y(v.x), z(v.y), w(b)   { }
    v4(float a, float b, v2 v) : x(a),   y(b),   z(v.x), w(v.y) { }
    v4(v3 v, float a)          : x(v.x), y(v.y), z(v.z), w(a)   { }
    v4(float a, v3 v)          : x(a),   y(v.x), z(v.y), w(v.z) { }
    
    v3 xyz()
    {
        return v3(x,y,z);
    }
};

struct v2i
{
    __m128i v;
    v2i() : v(_mm_set1_epi32(0)) {}
    v2i(int x, int y) : v(_mm_setr_epi32(x, y, 0, 0)) {}
    v2i(__m128i r) : v(r) {}
    int x() const { return _mm_extract_epi32(v, 0); }
    int y() const { return _mm_extract_epi32(v, 1); }
};

struct v3i
{
    __m128i v;
    v3i() : v(_mm_set1_epi32(0)) {}
    v3i(int x, int y, int z) : v(_mm_setr_epi32(x, y, z, 0)) {}
    v3i(__m128i r) : v(r) {}
    int x() const { return _mm_extract_epi32(v, 0); }
    int y() const { return _mm_extract_epi32(v, 1); }
    int z() const { return _mm_extract_epi32(v, 2); }
};

struct mat3
{
    float m[3][3];
    
    mat3() : m {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 0.0f} }
    {}
    
    mat3(float aa, float ab, float ac,
         float ba, float bb, float bc,
         float ca, float cb, float cc) :
    m { {aa, ab, ac},
        {ba, bb, bc},
        {ca, cb, cc} }
    {}
    
    const float *operator[](unsigned i) const
    {
        return m[i];
    }
    float *operator[](unsigned i)
    {
        return m[i];
    }
};

struct mat4
{
    float m[4][4];
    
    mat4() : m {
        {1.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f} }
    {}
    
    mat4(float aa, float ab, float ac, float ad,
         float ba, float bb, float bc, float bd,
         float ca, float cb, float cc, float cd,
         float da, float db, float dc, float dd) :
    m { {aa, ab, ac, ad},
        {ba, bb, bc, bd},
        {ca, cb, cc, cd},
        {da, db, dc, dd} }
    {}
    
    const float *operator[](unsigned i) const
    {
        return m[i];
    }
    float *operator[](unsigned i)
    {
        return m[i];
    }
};

template<u64 W, u64 H, u64 D>
struct BitCube
{
    static constexpr u64 N = W*H*D;
    static_assert(N % 8 == 0);
    u8 m_v[N/8] = {};
    
    void assign_bit(u64 n, bool v)
    {
        const u64 chunk = n >> 3;
        const u64 pos = n & 7;
        m_v[chunk] &= ~(1 << pos);
        m_v[chunk] |= (v << pos);
    }
    
    void set_bit(u64 n)
    {
        const u64 chunk = n >> 3;
        const u64 pos = n & 7;
        m_v[chunk] |= (1 << pos);
    }
    
    void clear_bit(u64 n)
    {
        const u64 chunk = n >> 3;
        const u64 pos = n & 7;
        m_v[chunk] &= ~(1 << pos);
    }
    
    bool is_bit_set(u64 x, u64 y, u64 z)
    {
        const u64 n = z*W*H + y*W + x;
        const u64 chunk = n >> 3;
        const u64 pos = n & 7;
        return m_v[chunk] & (1 << pos);
    }
    
    bool is_bit_set(u64 n)
    {
        const u64 chunk = n >> 3;
        const u64 pos = n & 7;
        return m_v[chunk] & (1 << pos);
    }
    
    bool all_bits_set()
    {
        // Assuming divisible by 8
        for(u32 i = 0; i < N/8; i++)
        {
            if(m_v[i] != 0xFF) return false;
        }
        return true;
    }
};

template<u64 N>
struct BitArray
{
    u64 v[N/64] = {};
    
    void set_bit(u64 n)
    {
        const u64 chunk = n / 64;
        const u64 pos = n & 63;
        v[chunk] |= (1ULL << pos);
    }

    void clear_bit(u64 n)
    {
        const u64 chunk = n / 64;
        const u64 pos = n & 63;
        u64 data = v[chunk];
        u64 clear_mask = 1;
        clear_mask = clear_mask << pos;
        clear_mask = ~clear_mask;
        data = data & clear_mask;
        v[chunk] = data;
    }
    
    bool is_bit_set(u64 n)
    {
        const u64 chunk = n / 64;
        const u64 pos = n & 63;
        return v[chunk] & (1ULL << pos);
    }
    
    bool all_bits_set()
    {
        static_assert(N % 64 == 0);
        u64 mask = u64(-1);
        for(u32 i = 0; i < N/64; i++)
        {
            mask &= v[i];
        }
        return !bool(~mask);
    }

    u64 tzcnt()
    {
        static_assert(N % 64 == 0);
        u64 result = 0;
        u64 i = 0;
        u64 tz;
        do
        {
            tz = _tzcnt_u64(v[i]);
            result += tz;
            i++;
        }
        while(tz == 64 && result < N);
        return result;
    }
    
    void clear()
    {
        memset(v, 0, sizeof(v));
    }
};

inline float squared(float a) { return a * a; }
inline u32 squared(u32 a) { return a * a; }
inline s32 squared(s32 a) { return a * a; }
float lerp(float a, float b, float t) { return ((1.0f - t) * a) + (t * b); }
float deg_to_rad(float a) { return a * (M_PI / 180.0f); }
float rad_to_deg(float a) { return a * (180.0f / M_PI); }
//float floor(float a) { return float(int(a)); }
//float ceil(float a) { return float(int(a - 1.0f)+1); }
//float ceil(float a) { return ceilf(a); }
inline void swap(s32& a, s32& b)
{
    s32 tmp = a;
    a = b;
    b = tmp;
}
inline void swap(u32& a, u32& b)
{
    u32 tmp = a;
    a = b;
    b = tmp;
}
inline void swap(f32& a, f32& b)
{
    f32 tmp = a;
    a = b;
    b = tmp;
}
inline bool is_power_of_2(s32 n)
{
    return (n & (n - 1)) == 0;
}
float inv_lerp(float a, float b, float t) { return (t - a) / (b - a); }
float remap(float a, float from_min, float from_max, float to_min, float to_max)
{
    float from_percent = inv_lerp(from_min, from_max, a);
    return lerp(to_min, to_max, from_percent);
}
int clamp(int a, int min, int max) { if(a < min) return min; if(a > max) return max; return a; }
/*
   float cos(float a) { return cosf(a); }
   float sin(float a) { return sinf(a); }
   float sqrt(float a) { return sqrtf(a); }
   float to_power(float base, float exponent) { return (float)pow(base, exponent); }
   int min(int a, int b) { return (a < b) ? a : b; }
   int min(int a, int b, int c) { return min(a, min(b, c)); }
   int max(int a, int b) { return (a > b) ? a : b; }
   int max(int a, int b, int c) { return max(a, max(b, c)); }
   float min(float a, float b) { return (a < b) ? a : b; }
   float min(float a, float b, float c) { return min(a, min(b, c)); }
   float max(float a, float b) { return (a > b) ? a : b; }
   float max(float a, float b, float c) { return max(a, max(b, c)); }
   float clamp(float a, float min, float max) { if(a < min) return min; if(a > max) return max; return a; }
   int floor(float a) { return (int)a; }
   float abs(float a) { return (a < 0.0f) ? -a : a; } 
   float average(float *values, int num)
   {
   float sum = 0.0f;
   for(int i = 0; i < num; i++) sum += values[i];
   return sum / (float)num;
   }
   float infinity()
   {
   return INFINITY;
   }
   */






static v2 operator+(v2 a, v2 b) { return v2(a.x + b.x, a.y + b.y); }
static v3 operator+(v3 a, v3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static v4 operator+(v4 a, v4 b) { return v4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }
static v2i operator+(v2i a, v2i b) { return _mm_add_epi32(a.v, b.v); }
static v3i operator+(v3i a, v3i b) { return _mm_add_epi32(a.v, b.v); }

static v2 operator-(v2 a, v2 b) { return v2(a.x - b.x, a.y - b.y); }
static v3 operator-(v3 a, v3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static v4 operator-(v4 a, v4 b) { return v4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }
static v2i operator-(v2i a, v2i b) { return _mm_sub_epi32(a.v, b.v); }
static v3i operator-(v3i a, v3i b) { return _mm_sub_epi32(a.v, b.v); }

static v2 operator-(v2 a) { return v2(-a.x, -a.y); }
static v3 operator-(v3 a) { return v3(-a.x, -a.y, -a.z); }
static v4 operator-(v4 a) { return v4(-a.x, -a.y, -a.z, -a.w); }
static v2i operator-(v2i a) { return _mm_sub_epi32(_mm_set1_epi32(0), a.v); }
static v3i operator-(v3i a) { return _mm_sub_epi32(_mm_set1_epi32(0), a.v); }

static bool operator==(v2 a, v2 b)
{
    return (a.x == b.x) && (a.y == b.y);
}

static bool operator==(v3i a, v3i b)
{
    assert(_mm_extract_epi32(a.v, 3) == 0);
    assert(_mm_extract_epi32(b.v, 3) == 0);
    __m128i m = _mm_xor_si128(a.v, b.v);
    m = _mm_cmpeq_epi32(m, _mm_set1_epi32(0));
    m = _mm_andnot_si128(m, _mm_set1_epi32(u32(-1)));
    return !bool(_mm_movemask_epi8(m));
}

static float dot(v2 a, v2 b) { return (a.x * b.x) + (a.y * b.y); }
static float dot(v3 a, v3 b) { return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }
static float dot(v4 a, v4 b) { return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w); }

static s32 dot(v3i a, v3i b)
{
    __m128i m = _mm_mullo_epi32(a.v, b.v);
    __m128i sum = _mm_add_epi32(
        m,
        _mm_add_epi32(
            _mm_shuffle_epi32(m, 0b00'11'10'01),
            _mm_add_epi32(
                _mm_shuffle_epi32(m, 0b01'00'11'10),
                _mm_shuffle_epi32(m, 0b10'01'00'11)
            )
        )
    );
    return _mm_cvtsi128_si32(sum);
}

static v2 operator*(v2 a, float scalar) { return v2(a.x * scalar, a.y * scalar); }
static v3 operator*(v3 a, float scalar) { return v3(a.x * scalar, a.y * scalar, a.z * scalar); }
static v4 operator*(v4 a, float scalar) { return v4(a.x * scalar, a.y * scalar, a.z * scalar, a.w * scalar); }

static v3i operator*(v3i a, s32 scalar)
{
    return _mm_mullo_epi32(a.v, _mm_set1_epi32(scalar));
}

static v2 operator*(float scalar, v2 a) { return v2(a.x * scalar, a.y * scalar); }
static v3 operator*(float scalar, v3 a) { return v3(a.x * scalar, a.y * scalar, a.z * scalar); }
static v4 operator*(float scalar, v4 a) { return v4(a.x * scalar, a.y * scalar, a.z * scalar, a.w * scalar); }

static v2 multiply_components(v2 a, v2 b) { return v2(a.x * b.x, a.y * b.y); }
static v3 multiply_components(v3 a, v3 b) { return v3(a.x * b.x, a.y * b.y, a.z * b.z); }
static v4 multiply_components(v4 a, v4 b) { return v4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w); }

static v2 operator/(v2 a, float scalar) { return v2(a.x / scalar, a.y / scalar); }
static v3 operator/(v3 a, float scalar) { return v3(a.x / scalar, a.y / scalar, a.z / scalar); }
static v4 operator/(v4 a, float scalar) { return v4(a.x / scalar, a.y / scalar, a.z / scalar, a.w / scalar); }

static v2 divide_components(v2 a, v2 b) { return v2(a.x / b.x, a.y / b.y); }
static v3 divide_components(v3 a, v3 b) { return v3(a.x / b.x, a.y / b.y, a.z / b.z); }
static v4 divide_components(v4 a, v4 b) { return v4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }

static v2 &operator+=(v2 &a, v2 b) { a.x += b.x; a.y += b.y; return a; }
static v3 &operator+=(v3 &a, v3 b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
static v4 &operator+=(v4 &a, v4 b) { a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w; return a; }

static v2 &operator-=(v2 &a, v2 b) { a.x -= b.x; a.y -= b.y; return a; }
static v3 &operator-=(v3 &a, v3 b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
static v4 &operator-=(v4 &a, v4 b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w; return a; }

static v2 &operator*=(v2 &a, float b) { a.x *= b; a.y *= b; return a; }
static v3 &operator*=(v3 &a, float b) { a.x *= b; a.y *= b; a.z *= b; return a; }
static v4 &operator*=(v4 &a, float b) { a.x *= b; a.y *= b; a.z *= b; a.w *= b; return a; }

static v2 &operator/=(v2 &a, float b) { a.x /= b; a.y /= b; return a; }
static v3 &operator/=(v3 &a, float b) { a.x /= b; a.y /= b; a.z /= b; return a; }
static v4 &operator/=(v4 &a, float b) { a.x /= b; a.y /= b; a.z /= b; a.w /= b; return a; }

static v2 reflect(v2 v, v2 n) { return v - 2.0f * dot(v, n) * n; }
static v3 reflect(v3 v, v3 n) { return v - 2.0f * dot(v, n) * n; }
static v4 reflect(v4 v, v4 n) { return v - 2.0f * dot(v, n) * n; }

static v2 lerp(v2 a, v2 b, float t) { return ((1.0f - t) * a) + (t * b); }
static v3 lerp(v3 a, v3 b, float t) { return ((1.0f - t) * a) + (t * b); }
static v4 lerp(v4 a, v4 b, float t) { return ((1.0f - t) * a) + (t * b); }

static int cube_idx(v3i a, int dim)
{
    // a.z*dim*dim + a.y*dim + a.x
    __m128i d = _mm_setr_epi32(1, dim, dim*dim, 0);
    __m128i r = _mm_mullo_epi32(a.v, d);
    __m128i r_shuf0 = _mm_shuffle_epi32(r, 0b00'11'10'01);
                                        __m128i r_shuf1 = _mm_shuffle_epi32(r, 0b01'00'11'10);
                                                                            r = _mm_add_epi32(r, r_shuf0);
                                                                            r = _mm_add_epi32(r, r_shuf1);
                                                                            return _mm_cvtsi128_si32(r);
}

static v2 lerp_components(v2 a, v2 b, v2 t)
{
    return multiply_components(v2(1.0f, 1.0f) - t, a) + multiply_components(t, b);
}

static v3 lerp_components(v3 a, v3 b, v3 t)
{
    return multiply_components(v3(1.0f, 1.0f, 1.0f) - t, a) + multiply_components(t, b);
}

static v4 lerp_components(v4 a, v4 b, v4 t)
{
    return multiply_components(v4(1.0f, 1.0f, 1.0f, 1.0f) - t, a) + multiply_components(t, b);
}

static v2 inv_lerp_components(v2 a, v2 b, v2 t) { return divide_components((t - a), (b - a)); }
static v3 inv_lerp_components(v3 a, v3 b, v3 t) { return divide_components((t - a), (b - a)); }
static v4 inv_lerp_components(v4 a, v4 b, v4 t) { return divide_components((t - a), (b - a)); }

static v2 remap(v2 a, v2 from_min, v2 from_max, v2 to_min, v2 to_max) { return lerp_components(to_min, to_max, inv_lerp_components(from_min, from_max, a)); }
static v3 remap(v3 a, v3 from_min, v3 from_max, v3 to_min, v3 to_max) { return lerp_components(to_min, to_max, inv_lerp_components(from_min, from_max, a)); }
static v4 remap(v4 a, v4 from_min, v4 from_max, v4 to_min, v4 to_max) { return lerp_components(to_min, to_max, inv_lerp_components(from_min, from_max, a)); }


// Gets the length of the vector
static float length(v2 v)
{
    return sqrtf(squared(v.x) + squared(v.y));
}
static float length(v3 v)
{
    return sqrtf(squared(v.x) + squared(v.y) + squared(v.z));
}
static float length(v4 v)
{
    return sqrtf(squared(v.x) + squared(v.y) + squared(v.z) + squared(v.w));
}

// Gets the squared length of this vector
static float length_squared(v2 v)
{
    return squared(v.x) + squared(v.y);
}
static float length_squared(v3 v)
{
    return squared(v.x) + squared(v.y) + squared(v.z);
}
static float length_squared(v4 v)
{
    return squared(v.x) + squared(v.y) + squared(v.z) + squared(v.w);
}

static v2 normalize(v2 v)
{
    float l = length(v);
    if(l == 0.0f) return v2();
    return v / l;
}
static v3 normalize(v3 v)
{
    float l = length(v);
    if(l == 0.0f) return v3();
    return v / l;
}
static v4 normalize(v4 v)
{
    float l = length(v);
    if(l == 0.0f) return v4();
    return v / l;
}

static v2 normalize_nonzero(v2 v) { return v / length(v); }
static v3 normalize_nonzero(v3 v) { return v / length(v); }
static v4 normalize_nonzero(v4 v) { return v / length(v); }

static float signed_distance_to_plane(v2 a, v2 p, v2 n)
{
    n = normalize(n);
    return dot(a, n) - dot(p, n);
}
static float signed_distance_to_plane(v3 a, v3 p, v3 n)
{
    n = normalize(n);
    return dot(a, n) - dot(p, n);
}

static float angle_between(v2 a, v2 b) { return (float)acos(dot(normalize(a), normalize(b))); }
static float angle_between(v3 a, v3 b) { return (float)acos(dot(normalize(a), normalize(b))); }
static float angle_between(v4 a, v4 b) { return (float)acos(dot(normalize(a), normalize(b))); }

// Returns this vector clamped by max length
static v2 clamp_length(v2 v, float max_length)
{
    float len = length(v);
    if(len > max_length)
    {
        v = v / len;
        v = v * max_length;
    }
    
    return v;
}
static v3 clamp_length(v3 v, float max_length)
{
    float len = length(v);
    if(len > max_length)
    {
        v = v / len;
        v = v * max_length;
    }
    
    return v;
}
static v4 clamp_length(v4 v, float max_length)
{
    float len = length(v);
    if(len > max_length)
    {
        v = v / len;
        v = v * max_length;
    }
    
    return v;
}

static v2 project_onto_plane(v2 a, v2 p, v2 n)
{
    float dist_from_plane = signed_distance_to_plane(a, p, n);
    return a - n * dist_from_plane;
}




// Returns a counter-clockwise normal to the given vector with the same length
static v2 find_ccw_normal(v2 a)
{
    return v2(-a.y, a.x);
}

static bool is_ccw(v2 a, v2 b, v2 c)
{
    float cross_p = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
    return cross_p > 0.0f;
}


// Returns this vector rotated by the angle in radians
static v2 rotate_vector(v2 a, float angle)
{
    v2 v;
    
    v.x = a.x * (float)cos(angle) - a.y * (float)sin(angle);
    v.y = a.x * (float)sin(angle) + a.y * (float)cos(angle);
    
    return v;
}

static v2 slerp(v2 a, v2 b, float t)
{
    float angle = angle_between(a, b);
    v2 t1 = (sinf((1.0f - t) * angle) / sinf(angle)) * a;
    v2 t2 = (sinf(t * angle) / sinf(angle)) * b;
    return t1 + t2;
}

// Returns the angle this vector is pointing at in radians. If the vector is
// pointing straight right, the angle is 0. If left, the angle is PI / 2.0f
// etc.
//            PI/2
//             |
//       PI <-- --> 0
//             |
//         -PI/2
static float angle(v2 a)
{
    return atan2f(a.y, a.x);
}



static v3 cross(v3 a, v3 b)
{
    v3 v;
    v.x = (a.y * b.z) - (a.z * b.y);
    v.y = (a.z * b.x) - (a.x * b.z);
    v.z = (a.x * b.y) - (a.y * b.x);
    return v;
}

static v3 hsv_to_rgb(v3 hsv)
{
    float hh, p, q, t, ff;
    int i;
    v3 rgb;
    
    if(hsv.s <= 0.0)
    {
        rgb.r = hsv.v;
        rgb.g = hsv.v;
        rgb.b = hsv.v;
        return rgb;
    }
    hh = hsv.h;
    if(hh >= 360.0f) hh = 0.0f;
    hh /= 60.0f;
    i = (int)hh;
    ff = hh - i;
    p = hsv.v * (1.0f - hsv.s);
    q = hsv.v * (1.0f - (hsv.s * ff));
    t = hsv.v * (1.0f - (hsv.s * (1.0f - ff)));
    
    switch(i)
    {
        case 0:
            rgb.r = hsv.v;
            rgb.g = t;
            rgb.b = p;
            break;
        case 1:
            rgb.r = q;
            rgb.g = hsv.v;
            rgb.b = p;
            break;
        case 2:
            rgb.r = p;
            rgb.g = hsv.v;
            rgb.b = t;
            break;

        case 3:
            rgb.r = p;
            rgb.g = q;
            rgb.b = hsv.v;
            break;
        case 4:
            rgb.r = t;
            rgb.g = p;
            rgb.b = hsv.v;
            break;
        case 5:
        default:
            rgb.r = hsv.v;
            rgb.g = p;
            rgb.b = q;
            break;
    }
    return rgb;     
}




// This function was made only for the matrix-vector multiplication
static float dot3v(const float *a, v3 b)
{
    return (a[0] * b.x) + (a[1] * b.y) + (a[2] * b.z);
}
static v3 operator*(const mat3 &lhs, v3 rhs)
{
    v3 result;
    result.x = dot3v(lhs[0], rhs);
    result.y = dot3v(lhs[1], rhs);
    result.z = dot3v(lhs[2], rhs);
    return result;
}
static mat3 operator*(const mat3 &lhs, const mat3 &rhs)
{
    mat3 product;
    for(unsigned row = 0; row < 3; row++)
    {
        for(unsigned col = 0; col < 3; col++)
        {
            float dot = 0.0f;
            for(unsigned i = 0; i < 3; i++)
            {
                dot += lhs.m[row][i] * rhs.m[i][col];
            }
            product[row][col] = dot;
        }
    }
    return product;
}
// This function was made only for the matrix-vector multiplication
static float dot4v(const float *a, v4 b)
{
    return (a[0] * b.x) + (a[1] * b.y) + (a[2] * b.z) + (a[3] * b.w);
}
static v4 operator*(const mat4 &lhs, v4 rhs)
{
    v4 result;
    
    result.x = dot4v(lhs[0], rhs);
    result.y = dot4v(lhs[1], rhs);
    result.z = dot4v(lhs[2], rhs);
    result.w = dot4v(lhs[3], rhs);
    
    return result;
}

static mat4 operator*(const mat4 &lhs, const mat4 &rhs)
{
    mat4 product;
    
    // Loop through each spot in the resulting matrix
    for(unsigned row = 0; row < 4; row++)
    {
        for(unsigned col = 0; col < 4; col++)
        {
            // Dot the row and column for the given slot
            float dot = 0.0f;
            for(unsigned i = 0; i < 4; i++)
            {
                dot += lhs.m[row][i] * rhs.m[i][col];
            }
            
            product[row][col] = dot;
        }
    }
    
    return product;
}




















// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][] 
// https://www.geeksforgeeks.org/adjoint-inverse-matrix/
static void get_cofactor(mat4 &a, mat4 &temp, int p, int q, int n) 
{ 
    int i = 0, j = 0; 
    
    // Looping for each element of the matrix 
    for(int row = 0; row < n; row++) 
    { 
        for(int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if(row != p && col != q) 
            { 
                temp[i][j++] = a[row][col]; 
                
                // Row is filled, so increase row index and 
                // reset col index 
                if(j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 

// Recursive function for finding determinant of matrix. 
// n is current dimension of a[][]
static float determinant(mat4 &A, int n) 
{ 
    if(n == 1) return A[0][0];
    
    float D = 0; // Initialize result 
    mat4 temp; // To store cofactors 
    float sign = 1;  // To store sign multiplier 
    
    // Iterate for each element of first row 
    for(int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        get_cofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 
        
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
    
    return D; 
} 

// Function to get adjoint of A[N][N] in adj[N][N]. 
static void adjoint(mat4 &A, mat4 &adj) 
{ 
    // temp is used to store cofactors of A[][] 
    int sign = 1;
    mat4 temp; 
    
    for(int i = 0; i < 4; i++) 
    { 
        for(int j = 0; j < 4; j++) 
        { 
            // Get cofactor of A[i][j] 
            get_cofactor(A, temp, i, j, 4); 
            
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i + j) % 2 == 0) ? 1 : -1; 
            
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = sign * determinant(temp, 4-1); 
        } 
    } 
} 

// Function to calculate and store inverse, returns false if 
// matrix is singular 
static mat4 inverse(mat4 A) 
{
    mat4 inverse;
    
    // Find determinant of A[][] 
    float det = determinant(A, 4); 
    if (det == 0) 
    { 
        return mat4(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); 
    } 
    
    // Find adjoint 
    mat4 adj; 
    adjoint(A, adj); 
    
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for(int i = 0; i < 4; i++) 
        for(int j = 0; j < 4; j++) 
        inverse[i][j] = adj[i][j] / det; 
    
    return inverse; 
}










static mat4 make_translation_matrix(v3 offset)
{
    mat4 result = 
    {
        1.0f, 0.0f, 0.0f, offset.x,
        0.0f, 1.0f, 0.0f, offset.y,
        0.0f, 0.0f, 1.0f, offset.z,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    return result;
}

static mat4 make_scale_matrix(v3 scale)
{
    mat4 result = 
    {
        scale.x, 0.0f, 0.0f, 0.0f,
        0.0f, scale.y, 0.0f, 0.0f,
        0.0f, 0.0f, scale.z, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    return result;
}

static mat4 make_x_axis_rotation_matrix(float radians)
{
    mat4 result = 
    {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, (float)cos(radians), (float)-sin(radians), 0.0f,
        0.0f, (float)sin(radians), (float)cos(radians), 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    return result;
}

static mat4 make_y_axis_rotation_matrix(float radians)
{
    mat4 result = 
    {
        (float)cos(radians), 0.0f, (float)sin(radians), 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        (float)-sin(radians), 0.0f, (float)cos(radians), 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    return result;
}

static mat4 make_z_axis_rotation_matrix(float radians)
{
    mat4 result = 
    {
        (float)cos(radians), (float)-sin(radians), 0.0f, 0.0f,
        (float)sin(radians), (float)cos(radians), 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    
    return result;
}

static mat3 make_axis_rotation_matrix(float rad, v3 u)
{
    float c = cosf(rad);
    float s = sinf(rad);
    mat3 result =
    {
        c + u.x*u.x*(1.0f - c),      u.x*u.y*(1.0f - c) - u.z*s,  u.x*u.z*(1.0f - c) + u.y*s,
        u.y*u.x*(1.0f - c) + u.z*s,  c + u.y*u.y*(1.0f - c),      u.y*u.z*(1.0f - c) - u.x*s,
        u.z*u.x*(1.0f - c) - u.y*s,  u.z*u.y*(1.0f - c) + u.x*s,  c + u.z*u.z*(1.0f - c)
    };
    return result;
}



static void seed_random(unsigned int n)
{
    srand(n);
}

static float random_01()
{
    return (float)rand() / (float)RAND_MAX;
}

static int random_range(int min, int max)
{
    int random_number = rand();
    
    int diff = max - min;
    int range = random_number % (diff + 1);
    int result = range + min;
    return result;
}

static float random_range(float min, float max)
{
    return lerp(min, max, random_01());
}

static v4 random_color()
{
    float hue = random_range(0.0f, 360.0f);
    float sat = 1.0f;
    float value = 1.0f;
    v3 hsv = { hue, sat, value };
    return v4( hsv_to_rgb(hsv), 1.0f);
}


static u32 random(u32 x)
{
    x ^= x << 13;
    x *= 182376581;
    x ^= x >> 17;
    x *= 783456103;
    x ^= x << 5;
    x *= 53523;
    return x;
}

static __m256i rand8(__m256i x)
{
    x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 13));
    x = _mm256_mullo_epi32(x, _mm256_set1_epi32(182376581));
    x = _mm256_xor_si256(x, _mm256_srli_epi32(x, 17));
    x = _mm256_mullo_epi32(x, _mm256_set1_epi32(783456103));
    x = _mm256_xor_si256(x, _mm256_slli_epi32(x, 5));
    x = _mm256_mullo_epi32(x, _mm256_set1_epi32(53523));
    return x;
}

static v2 random_gradient(s32 x, s32 y)
{
    u32 w = 8 * sizeof(unsigned);
    u32 s = w / 2;
    u32 a = x, b = y;
    a *= 3284157443;
    b ^= (a << s) | (a >> (w-s));
    b *= 1911520717;
    a ^= (b << s) | (b >> (w-s));
    a *= 2048419325;
    float random = a * (3.14159265f / ~(~0u >> 1));
    return v2(cosf(random), sinf(random));
}
static float smoothstep(float a0, float a1, float w)
{
    return (a1 - a0) * (3.0f - w * 2.0f) * w * w + a0;
}
static float perlin_noise(float x, float y)
{
    v2 p = v2(x, y);
    
    float x0 = x > 0.0f ? float(s32(x)) : float(s32(x - 1.0f));
    float x1 = x0 + 1;
    float y0 = y > 0.0f ? float(s32(y)) : float(s32(y - 1.0f));
    float y1 = y0 + 1;
    
    v2 d0 = p - v2(x0, y0);
    v2 d1 = p - v2(x1, y0);
    v2 d2 = p - v2(x0, y1);
    v2 d3 = p - v2(x1, y1);
    
    v2 g0 = random_gradient(s32(x0), s32(y0));
    v2 g1 = random_gradient(s32(x1), s32(y0));
    v2 g2 = random_gradient(s32(x0), s32(y1));
    v2 g3 = random_gradient(s32(x1), s32(y1));
    
    float i0 = smoothstep(dot(d0, g0), dot(d1, g1), x - x0);
    float i1 = smoothstep(dot(d2, g2), dot(d3, g3), x - x0);
    float i2 = smoothstep(i0, i1, y - y0);
    return i2;
}

f32 fract(f32 x)
{
    f32 f = std::floor(x);
    return x - f;
}

f32 approx_sin(f32 x)
{
    constexpr f32 pi = f32(M_PI);
    
    //f32 sign = x < 0 ? -1.0f : 1.0f;
    
    f32 sign = 1.0f;
    {
        f32 xx = x / (pi*2.0f);
        xx = fract(xx);
        sign = xx >= 0.5f ? -1.0f : 1.0f;
    }
    
    bool flip = false;
    {
        f32 xx = std::abs(x);
        xx = xx / pi;
        xx = fract(xx);
        flip = xx >= 0.5f;
    }
    
    x = std::abs(x);
    x = x / (pi*0.5f);
    x = fract(x);
    x *= pi*0.5f;
    
    if(flip)
    {
        x = pi*0.5f - x;
    }
    
    // [0, pi/2]  : x
    // [pi/2, pi] : 1.0f - x
    
    // [0, pi]    :  x
    // [pi, 2*pi] : -x
    
#define POW2(n) ((n)*(n))
#define POW3(n) ((n)*(n)*(n))
#define POW4(n) ((n)*(n)*(n)*(n))
    f32 eval_1 = x - POW3(x) / 6.0f;
    f32 eval_2 = (1.0f / 384.0f) * POW4(pi - 2.0f*x) - (1.0f/8.0f) * POW2(pi - 2.0f * x) + 1.0f;
    f32 eval = x < 0.6403 ? eval_1 : eval_2;
    return eval * sign;
#undef POW2
#undef POW3
#undef POW4
}

static inline __m256 approx_sin8(__m256 x)
{
    const __m256 pi       = _mm256_set1_ps(3.14159265f);
    const __m256 h_pi     = _mm256_set1_ps(0.5f * 3.14159265f);
    const __m256 inv_pi   = _mm256_set1_ps(1.0f / 3.14159265f);
    const __m256 inv_h_pi = _mm256_set1_ps(1.0f / (0.5f * 3.14159265f));
    const __m256 inv_2_pi = _mm256_set1_ps(1.0f / (2.0f * 3.14159265f));
    
    // Range reduce to [0, 2*pi] to check if we need to negate result.
    __m256 sign;
    {
        __m256 xx = _mm256_mul_ps(x, inv_2_pi);
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

static inline __m256 approx_cos8(__m256 x)
{
    const __m256 h_pi = _mm256_set1_ps(0.5f * 3.14159265f);
    return approx_sin8(_mm256_sub_ps(x, h_pi));
}

static inline __m256 lerp8(__m256 a, __m256 b, __m256 t)
{
    return _mm256_add_ps(
            _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(1.0f), t), a),
            _mm256_mul_ps(t, b)
            );
}

static inline __m256 pnoise8(__m256 x, __m256 y, __m256 z)
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

    const auto calc_dp = [x, y, z](const __m256 vx, const __m256 vy, const __m256 vz)
    {
#if 0
        // Alternate gradient function
        constexpr u32 num_sphere_points = 1 << 12;
        // Create noise
        __m256i i_noise = rand8(_mm256_castps_si256(_mm256_xor_ps(_mm256_xor_ps(vx, vy), vz)));
        i_noise = _mm256_and_si256(i_noise, _mm256_set1_epi32(num_sphere_points - 1));
        // Use random number to index points on a sphere.
        const __m256 noise = _mm256_cvtepi32_ps(i_noise);
        const __m256 u = _mm256_fmsub_ps(_mm256_set1_ps(2.0f / float(num_sphere_points - 1)), noise, _mm256_set1_ps(1.0f));
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
        const __m256 dx = _mm256_sub_ps(x, vx);
        const __m256 dy = _mm256_sub_ps(y, vy);
        const __m256 dz = _mm256_sub_ps(z, vz);
        const __m256 result = _mm256_fmadd_ps(gx, dx, _mm256_fmadd_ps(gy, dy, _mm256_mul_ps(dz, gz)));

#else
        __m256 result;
        {
            __m256i gxi = rand8(_mm256_castps_si256(_mm256_xor_ps(_mm256_xor_ps(vx, vy), vz)));
            gxi = _mm256_and_si256(gxi, _mm256_set1_epi32(0xFFFF));
            const __m256 gx = _mm256_mul_ps(_mm256_cvtepi32_ps(gxi), _mm256_set1_ps(1.0f/65535.0f));
            const __m256 dx = _mm256_sub_ps(x, vx);
            result = _mm256_mul_ps(dx, gx);
        }

        {
            __m256i gyi = rand8(
                    _mm256_mullo_epi32(
                        _mm256_castps_si256(_mm256_xor_ps(_mm256_xor_ps(vx, vy), vz)),
                        _mm256_set1_epi32(100)));
            gyi = _mm256_and_si256(gyi, _mm256_set1_epi32(0xFFFF));
            const __m256 gy = _mm256_mul_ps(_mm256_cvtepi32_ps(gyi), _mm256_set1_ps(1.0f/65535.0f));
            const __m256 dy = _mm256_sub_ps(y, vy);
            result = _mm256_fmadd_ps(dy, gy, result);
        }

        {
            const __m256 dz = _mm256_sub_ps(z, vz);
            __m256i gzi = rand8(
                    _mm256_mullo_epi32(
                        _mm256_castps_si256(_mm256_xor_ps(_mm256_xor_ps(vx, vy), vz)),
                        _mm256_set1_epi32(10000)));
            gzi = _mm256_and_si256(gzi, _mm256_set1_epi32(0xFFFF));
            const __m256 gz = _mm256_mul_ps(_mm256_cvtepi32_ps(gzi), _mm256_set1_ps(1.0f/65535.0f));
            result = _mm256_fmadd_ps(dz, gz, result);
        }
#endif

        return result;
    };

    // Smooth t.
    const __m256 dx0 = _mm256_sub_ps(x, x0);
    const __m256 dy0 = _mm256_sub_ps(y, y0);
    const __m256 dz0 = _mm256_sub_ps(z, z0);
    __m256 tx = _mm256_mul_ps(dx0, _mm256_mul_ps(dx0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dx0, _mm256_set1_ps(2.0f)))));
    __m256 ty = _mm256_mul_ps(dy0, _mm256_mul_ps(dy0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dy0, _mm256_set1_ps(2.0f)))));
    __m256 tz = _mm256_mul_ps(dz0, _mm256_mul_ps(dz0, _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(dz0, _mm256_set1_ps(2.0f)))));

    __m256 p0 = calc_dp(x0, y0, z0);
    __m256 p1 = calc_dp(x1, y0, z0);
    __m256 r0 = lerp8(p0, p1, tx);
    p0 = calc_dp(x0, y1, z0);
    p1 = calc_dp(x1, y1, z0);
    __m256 r1 = lerp8(p0, p1, tx);
    __m256 r2 = lerp8(r0, r1, ty);

    p0 = calc_dp(x0, y0, z1);
    p1 = calc_dp(x1, y0, z1);
    r0 = lerp8(p0, p1, tx);
    p0 = calc_dp(x0, y1, z1);
    p1 = calc_dp(x1, y1, z1);
    r1 = lerp8(p0, p1, tx);
    __m256 r3 = lerp8(r0, r1, ty);
    
    __m256 r4 = lerp8(r2, r3, tz);

    return r4;
}

#pragma warning(default:4201)
#pragma warning(default:4505)

