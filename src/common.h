
#pragma once

#include "win32_crt.h"


#if defined(__clang__)
    #define MAYBE_UNUSED __attribute__((unused))
#elif defined(__GNUC__) || defined(__GNUG__)
    #define MAYBE_UNUSED __attribute__((unused))
#elif defined(_MSC_VER)
    #define MAYBE_UNUSED 
#endif


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
typedef const int s32_m;
typedef long long s64_m;
typedef float f32_m;
typedef double f64_m;

#define ARRAY_COUNT(N) (sizeof(N) / sizeof((N)[0]))

#define KB(N) (N * 1024ULL)
#define MB(N) (N * 1024ULL * 1024ULL)
#define GB(N) (N * 1024ULL * 1024ULL * 1024ULL)

