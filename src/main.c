
#include "common.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

INT WINAPI WinMainCRTStartup()
{
    HMODULE instance = GetModuleHandle(0);
    (void)instance;

    int big[4096] = {0};
    big[0] = 42;

    memset(big, 0, sizeof(big));

    volatile u64 s = 1;
    volatile u64 u = 1;

    s += s;
    s -= s;
    s *= s;
    s /= s;
    s %= s;
    s >>= 33;
    s <<= 33;

    u += u;
    u -= u;
    u *= u;
    u /= u;
    u %= u;
    u >>= 33;
    u <<= 33;

    float f = 1000.0f;
    double d = 1000000000.0;

    s32 i32f = (s32)f;
    s32 i32d = (s32)d;
    s32 u32f = (s32)f;
    s32 u32d = (s32)d;

    s64 i64f = (s64)f;
    s64 i64d = (s64)d;
    u64 u64f = (u64)f;
    u64 u64d = (u64)d;

    f = (float)i32f;
    d = (double)i32d;
    f = (float)u32f;
    d = (double)u32d;

    f = (float)i64f;
    d = (double)i64d;
    f = (float)u64f;
    d = (double)u64d;

    return 0;
}

