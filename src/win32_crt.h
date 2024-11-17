
#pragma once

// Avoid C runtime library
// https://hero.handmade.network/forums/code-discussion/t/94-guide_-_how_to_avoid_c_c++_runtime_on_windows
extern int _fltused;

void *memset(void *dst, int c, size_t count);
void *memcpy(void *dst, const void *src, size_t count);

