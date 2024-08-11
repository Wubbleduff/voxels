
#include "common.h"
#include "math.h"
#include "platform.h"
#include "game_state.h"

#include "graphics.h"

INTERNAL inline void make_world_to_camera_space_mtx(mtx4x4* result, const v3 pos, f32 pitch_turns, f32 yaw_turns)
{
    mtx4x4 y_rot_mtx;
    make_y_axis_rotation_mtx(&y_rot_mtx, -yaw_turns);

    mtx4x4 x_rot_mtx;
    make_x_axis_rotation_mtx(&x_rot_mtx, -pitch_turns);

    mtx4x4 translation_mtx;
    make_translation_mtx(&translation_mtx, v3_scale(pos, -1.0f));

    mtx4x4 world_to_cam_mtx_temp;
    mtx4x4_mul(&world_to_cam_mtx_temp, &y_rot_mtx, &translation_mtx);
    mtx4x4_mul(result, &x_rot_mtx, &world_to_cam_mtx_temp);
}

void do_one_frame(
        struct GameState* next_game_state,
        struct GameState* prev_game_state,
        struct MemoryArena* memory_arena)
{
    (void)memory_arena;

    v3 player_pos = prev_game_state->player_pos;
    f32_m player_pitch_turns = prev_game_state->player_pitch_turns;
    f32_m player_yaw_turns = prev_game_state->player_yaw_turns;

    {

        // Camera controls
        player_pitch_turns += (f32)is_key_down(KB_I) * 0.001f;
        player_pitch_turns -= (f32)is_key_down(KB_K) * 0.001f;
        player_yaw_turns += (f32)is_key_down(KB_J) * 0.001f;
        player_yaw_turns -= (f32)is_key_down(KB_L) * 0.001f;

        if(is_fps_mode())
        {
            s32_m mouse_screen_dx;
            s32_m mouse_screen_dy;
            get_mouse_delta(&mouse_screen_dx, &mouse_screen_dy);
            player_pitch_turns -= (f32)(mouse_screen_dy) * 0.0005f;
            player_yaw_turns -= (f32)(mouse_screen_dx) * 0.0005f;
        }

        mtx4x4 world_to_cam_mtx;
        make_world_to_camera_space_mtx(&world_to_cam_mtx, player_pos, player_pitch_turns, player_yaw_turns);

        v3 I_cam = {.m={world_to_cam_mtx.m[0*4 + 0], world_to_cam_mtx.m[0*4 + 1], world_to_cam_mtx.m[0*4 + 2]}};
        v3 J_cam = {.m={world_to_cam_mtx.m[1*4 + 0], world_to_cam_mtx.m[1*4 + 1], world_to_cam_mtx.m[1*4 + 2]}};
        v3 K_cam = {.m={world_to_cam_mtx.m[2*4 + 0], world_to_cam_mtx.m[2*4 + 1], world_to_cam_mtx.m[2*4 + 2]}};

        f32 speed = 1.0f;
        player_pos = v3_add(player_pos, v3_scale(K_cam, -speed * (f32)is_key_down(KB_W)));
        player_pos = v3_add(player_pos, v3_scale(K_cam,  speed * (f32)is_key_down(KB_S)));

        player_pos = v3_add(player_pos, v3_scale(I_cam,  speed * (f32)is_key_down(KB_D)));
        player_pos = v3_add(player_pos, v3_scale(I_cam, -speed * (f32)is_key_down(KB_A)));
        
        player_pos = v3_add(player_pos, v3_scale(J_cam,  speed * (f32)is_key_down(KB_SPACE)));
        player_pos = v3_add(player_pos, v3_scale(J_cam, -speed * (f32)is_key_down(KB_LCTRL)));

        next_game_state->player_pos = player_pos;
        next_game_state->player_pitch_turns = player_pitch_turns;
        next_game_state->player_yaw_turns = player_yaw_turns;
    }

    for(u64_m i = 0; i < 10; i++)
    {
        if(prev_game_state->terrain_loaded[i])
        {
            struct Terrain* src_terrain = prev_game_state->terrain_lod_list + i;
            struct Terrain* dst_terrain = next_game_state->terrain_lod_list + i;
            copy_terrain(dst_terrain, src_terrain);
            next_game_state->terrain_loaded[i] = 1;
        }

        if(!next_game_state->terrain_loaded[i])
        {
            struct Terrain* dst_terrain = next_game_state->terrain_lod_list + i;
            generate_terrain(
                &dst_terrain->num_vertices,
                dst_terrain->vx,
                dst_terrain->vy,
                dst_terrain->vz,
                dst_terrain->nx,
                dst_terrain->ny,
                dst_terrain->nz,
                &dst_terrain->num_indices,
                dst_terrain->indices,
                player_pos.x,
                player_pos.y,
                player_pos.z,
                (u8)i
            );
            next_game_state->terrain_loaded[i] = 1;
        }
    }
    next_game_state->terrain_active_lod = prev_game_state->terrain_active_lod;
    if(is_key_toggled_down(KB_N))
    {
        next_game_state->terrain_active_lod++;
        next_game_state->terrain_active_lod %= 10;
    }
    if(is_key_toggled_down(KB_P))
    {
        next_game_state->terrain_active_lod = next_game_state->terrain_active_lod == 0 ? 9 : next_game_state->terrain_active_lod - 1;
    }
}

void draw_game_state(struct GameState* game_state, struct MemoryArena* memory_arena)
{
    /*
     * Derivation for 3D perspective projection matrix
     *        
     *                                                    /
     *                                                  /
     *                                                /
     *                                              /
     *                                            /
     *                                          /
     *                                        /
     *                                      /
     *                                    /
     *                                  /
     *                                /
     *                              /
     *                            /
     *                          /
     *                        /
     *                      / |               V
     *                    /   |             .>+---------------+
     *                  /     |       -----/ D|               |
     *         Y      /       |  ----/        |               |
     *         ^    /       --R-/             |               |
     *         |  /    ----/  |               |               |
     *         |/  ---/       |               |               |
     *  Z <----C../           |----> N        |               |
     *  F     / \             |               |               |
     *       /    \           |               +---------------+
     *      V       \         |
     *     X          \       |
     *                  \     |
     *                    \   |
     *                      \ |
     *                        \
     *                          \
     *                            \
     *                              \
     *                                \
     *                                  \
     *                                    \
     *                                      \
     *                                        \
     *                                          \
     *                                            \
     *                                              \
     *                                                \
     *                                                  \
     *                                                    \
     *        
     *        
     *        
     * C : 3D camera point (Assume the camera is at the origin - (0, 0)
     * F : Normalized camera forward vector (The camera looks along the -Z axis. For something to be seen, it must be more -Z than the near plane.)
     *     In camera space, this will be (0, 0, 1)
     * n : Camera's near plane distance. For something to be seen, it must have a Z coordinate < -n.
     * V : Vertex to be projected
     *
     * The goal is to intersect the ray from the origin to the vertex with the near plane.
     *
     * Q : Point along the ray (solve for intersection)
     * Q = C + unit(V)*t, but since C is just the origin
     * Q = unit(V)*t
     *
     * Find the plane equation:
     * N : The near plane normal (In camera space, this will be (0, 0, -1)
     * S : a point on the near plane.
     * S = N*n
     * 
     * New plane equation : (P - S) * N = 0
     * We want to find where a point on the ray is equal to 0, so plug in Q for P:
     * (Q - S) * N = 0
     *
     * Expand:
     *
     * (unit(V)*t - S) * N = 0
     *
     * Solve for t:
     *
     * unit(V)*t*N - S*N = 0
     * t = S*N / (unit(V)*N)
     *
     * Plug t back in to the ray equation:
     *
     * Q = unit(V)*t
     * Q = unit(V)*(S*N / (unit(V)*N))
     *
     * Find in terms of V
     *
     * Q = unit(V)*(S*N / (unit(V)*N))
     *
     * Q = unit(V)*S*N
     *     -----------
     *     (unit(V)*N)
     *
     * Q = unit(V)*(N*n)*N
     *     ---------------
     *       (unit(V)*N)
     *      
     * Q = N*n*N*unit(V)
     *     -------------
     *      (unit(V)*N)
     *
     * Q = N*n*N*(V / ||V*V||)
     *     -------------------
     *      ((V / ||V*V||)*N)
     *
     * Q = N*n*N*V
     *     -------
     *     (V*N)
     *
     * Q = N*N*n*V
     *     -------
     *     (V*N)
     *
     * In 3D:
     * Q = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
     *     -----------------------------------  *  V
     *       (V_x*N_x + V_y*N_y + V_z*N_z)
     *
     * Q_x = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
     *       -----------------------------------  *  V_x
     *         (V_x*N_x + V_y*N_y + V_z*N_z)
     *
     * Q_y = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
     *       -----------------------------------  *  V_y
     *         (V_x*N_x + V_y*N_y + V_z*N_z)
     *
     * Q_z = (N_x*N_x*n + N_y*N_y*n * N_z*N_z*n)
     *       -----------------------------------  *  V_z
     *         (V_x*N_x + V_y*N_y + V_z*N_z)
     * 
     * Assume our object has been translated to camera space. In this case, N = (0, 0, -1) (Right-handed coordinate system)
     * N_x = 0
     * N_y = 0
     * N_z = -1
     *
     * Q_x =   n
     *       ------ * V_x
     *       (-V_z)
     *
     * Q_y =   n
     *       ------ * V_y
     *       (-V_z)
     *
     * Q_z =   n
     *       ------ * V_z
     *       (-V_z)
     *
     *
     * So now we have the point Q in camera space where Q is V perspective projected onto the near plane.
     * Our goal is to find Q_p in NDC space. So, we need to divide X and Y by the camera width and height.
     *
     * C_w : camera width
     * C_h : camera height
     *
     * Q_px = Q_x / C_w = n / C_w
     *                    ------- * V_x
     *                    (-V_z)        
     *
     * Q_py = Q_y / C_h = n / C_h
     *                    ------- * V_y
     *                    (-V_z)
     *
     * Q_pz should be between -1 and 1, so we need to remap it from [-n/-V_z, -f/-V_z] to [-1, 1]
     * f : far plane dist
     *
     * After dividing Q_pz by -V_z, we want the result to be -1 when V_z = -n and 1 when V_z = -f
     * Q_pz / -V_z = -1, V_z = -n
     * Q_pz / -V_z =  1, V_z = -f
     * What are the V_z values _before_ dividing?
     * Z = -n when V is on the near plane
     * Z =  f when V is on the far plane
     * For the coefficients, they must perform a remap operation from [-n, -f] to [-n,  f]
     *      z - (-n)
     * t = -----------
     *     (-f) - (-n)
     * Result: -n + t * (f - (-n))
     * ->
     *       z - (-n)
     * -n + ----------- * (f + n)
     *      (-f) - (-n)
     * ->
     *      z + n
     * -n + ----- * (f + n)
     *      n - f
     * ->
     * -n * (n - f)   (z + n)
     * ------------ + ------- * (f + n)
     *     n - f       n - f
     * -> 
     * -n * (n - f) + (z + n) * (f + n)
     * ---------------------------------
     *             n - f    
     * -> 
     * -n^2 + n*f  +  z*f + z*n + n*f + n^2
     * ---------------------------------
     *             n - f    
     * -> 
     * n*f  +  z(f + n) + n*f
     * -------------------------
     *          n - f    
     * -> 
     * (f + n) * z + 2*n*f
     * --------------------
     *         n - f
     * -> 
     * (f + n)        2 * n * f
     * -------- * z + ----------
     *  n - f           n - f
     * 
     *
     * https://www.desmos.com/calculator/krbxjns133
     *       
     * Now, define as a matrix (keep in mind we will be dividing by the W component after matrix multiplication):
     *
     * |      |   |  n                                  |   |      |
     * | Q_px |   | ---      0         0          0     |   | V_px |
     * |      |   | C_w                                 |   |      |
     * |      |   |                                     |   |      |
     * |      |   |          n                          |   |      |
     * | Q_py |   |  0      ---        0          0     |   | V_py |
     * |      | = |         C_h                         | * |      |
     * |      |   |                                     |   |      |
     * |      |   |                  f + n      2*n*f   |   |      |
     * | Q_pz |   |  0       0      -------    -------  |   | V_pz |
     * |      |   |                  n - f      n - f   |   |      |
     * |      |   |                                     |   |      |
     * | Q_pw |   |  0       0        -1          0     |   | V_pw |
     *
     * ->
     *
     * | Q_px |   |          (n / C_w) * V_x                     |
     * | Q_py | = |          (n / C_h) * V_y                     |
     * | Q_pz |   |  (f + n) / (n - f) * V_z + (2*n*f) / (n - f) |
     * | Q_pw |   |                 -1 * V_z                     | <-- Will be dividing all the terms by -V_z
     *
     */

    mtx4x4 world_to_cam_mtx;
    make_world_to_camera_space_mtx(&world_to_cam_mtx, game_state->player_pos, game_state->player_pitch_turns, game_state->player_yaw_turns);

    f32 aspect_ratio = get_screen_aspect_ratio();

    f32 n = 0.1f;
    f32 f = 1000.0f;
    f32 C_w = 0.125f;
    f32 C_h = C_w * aspect_ratio;
    mtx4x4 proj_mtx = {
        .m = {
            //    X         Y             Z                          W
            n / C_w,     0.0f,           0.0f,                      0.0f,
               0.0f,  n / C_h,           0.0f,                      0.0f,
               0.0f,     0.0f,    (f + n) / (n - f),  (2.0f * n * f) / (n - f),
               0.0f,     0.0f,          -1.0f,                      0.0f,
        }
    };

    mtx4x4 camera_and_clip_mtx;
    mtx4x4_mul(&camera_and_clip_mtx, &proj_mtx, &world_to_cam_mtx);

    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));
    struct Terrain* terrain = game_state->terrain_lod_list + game_state->terrain_active_lod;
    terrain_mesh->num_vertices = terrain->num_vertices;
    ARRAY_COPY(terrain_mesh->vx, terrain->vx, terrain->num_vertices);
    ARRAY_COPY(terrain_mesh->vy, terrain->vy, terrain->num_vertices);
    ARRAY_COPY(terrain_mesh->vz, terrain->vz, terrain->num_vertices);
    ARRAY_COPY(terrain_mesh->nx, terrain->nx, terrain->num_vertices);
    ARRAY_COPY(terrain_mesh->ny, terrain->ny, terrain->num_vertices);
    ARRAY_COPY(terrain_mesh->nz, terrain->nz, terrain->num_vertices);
    terrain_mesh->num_indices = terrain->num_indices;
    ARRAY_COPY(terrain_mesh->indices, terrain->indices, terrain->num_indices);
    f32 pos[3] = {0};
    draw_Mesh_1M(terrain_mesh, &camera_and_clip_mtx, 1, pos + 0, pos + 1, pos + 2);

    // Draw world basis.
    {
        v3 s = {.m={0.0f, 0.0f, 0.0f}};
        v3 e = {.m={1.0f, 0.0f, 0.0f}};
        debug_draw_line(&camera_and_clip_mtx, s, e, e);
        e.x = 0.0f;
        e.y = 1.0f;
        e.z = 0.0f;
        debug_draw_line(&camera_and_clip_mtx, s, e, e);
        e.x = 0.0f;
        e.y = 0.0f;
        e.z = 1.0f;
        debug_draw_line(&camera_and_clip_mtx, s, e, e);
    }
}


