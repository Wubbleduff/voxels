
#pragma once

#include "common.h"
#include "math.h"
#include "terrain.h"

// POD game state structure intended to hold all state for a particular frame of the game. Should be able to recreate a point in the game from 1 instance of this.
struct GameState
{
    v3 player_pos;
    f32 player_pitch_turns;
    f32 player_yaw_turns;

    struct Terrain terrain;
    struct TerrainProgress terrain_progress;
};

INTERNAL void reset_game_state(struct GameState* game_state)
{
    game_state->player_pos = v3_zero();
    game_state->player_pitch_turns = 0.0f;
    game_state->player_yaw_turns = 0.0f;

    reset_terrain(&game_state->terrain);
    memset(&game_state->terrain_progress, 0, sizeof(struct TerrainProgress));
}

