
#pragma once

#include "common.h"






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Memory

struct MemoryArena;
struct MemoryArena memory_arena_init(void* base, const u64 cap);
void memory_arena_reset(struct MemoryArena* arena);
void* memory_arena_allocate(const char* file, const s32 line, struct MemoryArena* arena, const u64 req_bytes);
void* memory_arena_allocate_zeroed(const char* file, const s32 line, struct MemoryArena* arena, const u64 req_bytes);
#define MEMORY_ARENA_ALLOCATE(arena, req_bytes) memory_arena_allocate(__FILE__, __LINE__, (arena), (req_bytes))
#define MEMORY_ARENA_ALLOCATE_ZEROED(arena, req_bytes) memory_arena_allocate_zeroed(__FILE__, __LINE__, (arena), (req_bytes))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Input

enum KeyboardKey
{
    KB_NOT_SUPPORTED,
    KB_ESCAPE,
    KB_SPACE,
    KB_LCTRL,
    KB_0,
    KB_1,
    KB_2,
    KB_3,
    KB_4,
    KB_5,
    KB_6,
    KB_7,
    KB_8,
    KB_9,
    KB_A,
    KB_B,
    KB_C,
    KB_D,
    KB_E,
    KB_F,
    KB_G,
    KB_H,
    KB_I,
    KB_J,
    KB_K,
    KB_L,
    KB_M,
    KB_N,
    KB_O,
    KB_P,
    KB_Q,
    KB_R,
    KB_S,
    KB_T,
    KB_U,
    KB_V,
    KB_W,
    KB_X,
    KB_Y,
    KB_Z,

    NUM_KEYBOARD_KEYS,
};
u32 is_key_down(enum KeyboardKey k);
u32 is_key_toggled_down(enum KeyboardKey k);
u32 is_fps_mode();
void get_mouse_delta(s32* x, s32* y);
void show_cursor(u32 show);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
