
#pragma once

#include "win32_crt.h"

#include <immintrin.h>


#if defined(__clang__)
    #define MAYBE_UNUSED __attribute__((unused))
#elif defined(__GNUC__) || defined(__GNUG__)
    #define MAYBE_UNUSED __attribute__((unused))
#elif defined(_MSC_VER)
    #define MAYBE_UNUSED 
#endif

#ifdef INTERNAL
#error "INTERNAL" already defined.
#endif
#define INTERNAL MAYBE_UNUSED static

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;
typedef char s8;
typedef short s16;
typedef int s32;
typedef long long s64;
typedef float f32;
typedef double f64;

_Static_assert(sizeof(u8) == 1, "Unexpected type size.");
_Static_assert(sizeof(u16) == 2, "Unexpected type size.");
_Static_assert(sizeof(u32) == 4, "Unexpected type size.");
_Static_assert(sizeof(u64) == 8, "Unexpected type size.");
_Static_assert(sizeof(s8) == 1, "Unexpected type size.");
_Static_assert(sizeof(s16) == 2, "Unexpected type size.");
_Static_assert(sizeof(s32) == 4, "Unexpected type size.");
_Static_assert(sizeof(s64) == 8, "Unexpected type size.");
_Static_assert(sizeof(f32) == 4, "Unexpected type size.");
_Static_assert(sizeof(f64) == 8, "Unexpected type size.");

#define u8_MAX 0xFF
#define u16_MAX 0xFFFF
#define u32_MAX 0xFFFFFFFF
#define u64_MAX 0xFFFFFFFFFFFFFFFF
#define s8_MAX 0x7F
#define s16_MAX 0x7FFF
#define s32_MAX 0x7FFFFFFF
#define s64_MAX 0x7FFFFFFFFFFFFFFF

void assert_fn(const char* file, int line, const u64 c, const char* msg);
#define ASSERT(c, msg) assert_fn(__FILE__, __LINE__, (const u64)(c), (msg))

#define ARRAY_COUNT(N) (sizeof(N) / sizeof((N)[0]))

#define COPY(DST, SRC, NUM) \
    do { \
        _Static_assert(sizeof((DST)[0]) == sizeof((SRC)[0]), "Array element sizes do not match."); \
        ASSERT(NUM <= ARRAY_COUNT(SRC), "Number of elements is larger than array size."); \
        ASSERT(NUM <= ARRAY_COUNT(DST), "Number of elements is larger than array size."); \
        memcpy((DST), (SRC), sizeof((SRC)[0]) * NUM); \
    } while(0)

#define KB(N) (N * 1024ULL)
#define MB(N) (N * 1024ULL * 1024ULL)
#define GB(N) (N * 1024ULL * 1024ULL * 1024ULL)


static inline u64 cstr_len(const char* s)
{
    u64 n = 0;
    while(*s++) n++;
    return n;
}

typedef struct
{
    // Last byte should always be 0 so we can trivially cast to C-string.
    _Alignas(32) char buf[32];
} StringBuf32;

static inline const char* StringBuf32_to_cstr(const StringBuf32* s)
{
    return s->buf;
}

static inline StringBuf32 u32_to_StringBuf32(u32 n)
{
    StringBuf32 result;
    memset(result.buf, 0, sizeof(result));
    u64 size = 0;
    while(n && size < 31)
    {
        const u32 digit = n % 10;
        result.buf[size++] = (char)('0' + digit);
        n /= 10;
    }
    for(u64 i = 0; i < size / 2; i++)
    {
        char tmp = result.buf[i];
        result.buf[i] = result.buf[size - i - 1];
        result.buf[size - i - 1] = tmp;
    }
    return result;
}

