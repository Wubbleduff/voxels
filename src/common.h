
#pragma once

#include "win32_crt.h"


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


typedef const unsigned char u8;
typedef const unsigned short u16;
typedef const unsigned int u32;
typedef const unsigned long long u64;
typedef const char s8;
typedef const short s16;
typedef const int s32;
typedef const long long s64;
typedef const float f32;
typedef const double f64;

typedef unsigned char u8_m;
typedef unsigned short u16_m;
typedef unsigned int u32_m;
typedef unsigned long long u64_m;
typedef char s8_m;
typedef short s16_m;
typedef int s32_m;
typedef long long s64_m;
typedef float f32_m;
typedef double f64_m;

#define MAX_U32 ((u32)(-1))

#define ARRAY_COUNT(N) (sizeof(N) / sizeof((N)[0]))

#define KB(N) (N * 1024ULL)
#define MB(N) (N * 1024ULL * 1024ULL)
#define GB(N) (N * 1024ULL * 1024ULL * 1024ULL)


static inline u64 cstr_len(const char* s)
{
    u64_m n = 0;
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

