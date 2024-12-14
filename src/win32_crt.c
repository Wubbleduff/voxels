
#include "common.h"

#define WIN32_LEAN_AND_MEAN 
#include <windows.h>

// Avoid C runtime library
// https://hero.handmade.network/forums/code-discussion/t/94-guide_-_how_to_avoid_c_c++_runtime_on_windows

int _fltused = 0;

#pragma function(memset)
void *memset(void *dst, int c, size_t count)
{
    char *bytes = (char *)dst;
    while (count--)
    {
        *bytes++ = (char)c;
    }
    return dst;
}

#pragma function(memcpy)
void *memcpy(void *dst, const void *src, size_t count)
{
    //char *dst8 = (char *)dst;
    //const char *src8 = (const char *)src;
    //while (count--)
    //{
    //    *dst8++ = *src8++;
    //}
    //return dst;

    if(count < 32)
    {
        for(u64 i = 0; i < count; i++)
        {
            ((u8*)dst)[i] = ((const u8*)src)[i];
        }
    }
    else
    {
        for(u64 i = 0; i < count; i += 32)
        {
            _mm256_storeu_si256((__m256i*)((u8*)dst + i), _mm256_loadu_si256((__m256i*)((const u8*)src + i)));
        }
    }

    for(u64 i = count & ~31; i < count; i++)
    {
        ((u8*)dst)[i] = ((const u8*)src)[i];
    }
    return (u8*)dst + count;
}

