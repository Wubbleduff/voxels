
#pragma once

#include "common.h"
#include "math.h"
#include "terrain.h"

// POD game state structure intended to hold all state for a particular frame of the game. Should be able to recreate a point in the game from 1 instance of this.
struct GameState
{
    v3 player_pos;
    f32_m player_pitch_turns;
    f32_m player_yaw_turns;

    u32_m terrain_generated;
    struct Terrain terrain;
};

INTERNAL void reset_game_state(struct GameState* game_state)
{
    game_state->player_pos = v3_zero();
    game_state->player_pitch_turns = 0.0f;
    game_state->player_yaw_turns = 0.0f;

    game_state->terrain_generated = 0;
    reset_terrain(&game_state->terrain);
}

