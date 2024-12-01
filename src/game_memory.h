
#include "game_state.h"
#include "platform_win32_opengl.h"

#pragma once

struct GameMemory
{
    struct Platform_Win32 platform_state;
    struct GameState game_state_A;
    struct GameState game_state_B;
};

