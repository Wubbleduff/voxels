
#include "common.h"
#include "platform.h"
#include "graphics.h"
#include "math.h"

#include "terrain.h"

void do_one_frame(struct MemoryArena* memory_arena)
{
    static v3 cam_pos = {.m = {0.0f, 1.0f, 3.0f}};

    struct Mesh_1M* terrain_mesh = (struct Mesh_1M*)MEMORY_ARENA_ALLOCATE(memory_arena, sizeof(struct Mesh_1M));

    generate_terrain(
        &terrain_mesh->num_vertices,
        terrain_mesh->vx,
        terrain_mesh->vy,
        terrain_mesh->vz,
        terrain_mesh->nx,
        terrain_mesh->ny,
        terrain_mesh->nz,
        &terrain_mesh->num_indices,
        terrain_mesh->indices,
        cam_pos.x,
        cam_pos.y,
        cam_pos.z
    );


    {
        // Camera controls
        static f32_m cam_pitch_turns = 0.0f;
        cam_pitch_turns += (f32)is_key_down(KB_I) * 0.001f;
        cam_pitch_turns -= (f32)is_key_down(KB_K) * 0.001f;


        static f32_m cam_yaw_turns = 0.0f;
        cam_yaw_turns += (f32)is_key_down(KB_J) * 0.001f;
        cam_yaw_turns -= (f32)is_key_down(KB_L) * 0.001f;

        if(is_fps_mode())
        {
            s32_m mouse_screen_dx;
            s32_m mouse_screen_dy;
            get_mouse_delta(&mouse_screen_dx, &mouse_screen_dy);
            cam_pitch_turns -= (f32)(mouse_screen_dy) * 0.0005f;
            cam_yaw_turns -= (f32)(mouse_screen_dx) * 0.0005f;
        }


        mtx4x4 y_rot_mtx;
        make_y_axis_rotation_mtx(&y_rot_mtx, -cam_yaw_turns);

        mtx4x4 x_rot_mtx;
        make_x_axis_rotation_mtx(&x_rot_mtx, -cam_pitch_turns);

        mtx4x4 translation_mtx;
        make_translation_mtx(&translation_mtx, v3_scale(cam_pos, -1.0f));

        mtx4x4 world_to_cam_mtx_temp;
        mtx4x4_mul(&world_to_cam_mtx_temp, &y_rot_mtx, &translation_mtx);

        mtx4x4 world_to_cam_mtx;
        mtx4x4_mul(&world_to_cam_mtx, &x_rot_mtx, &world_to_cam_mtx_temp);

        v3 I_cam = {.m={world_to_cam_mtx.m[0*4 + 0], world_to_cam_mtx.m[0*4 + 1], world_to_cam_mtx.m[0*4 + 2]}};
        v3 J_cam = {.m={world_to_cam_mtx.m[1*4 + 0], world_to_cam_mtx.m[1*4 + 1], world_to_cam_mtx.m[1*4 + 2]}};
        v3 K_cam = {.m={world_to_cam_mtx.m[2*4 + 0], world_to_cam_mtx.m[2*4 + 1], world_to_cam_mtx.m[2*4 + 2]}};

        f32 speed = 1.0f;
        cam_pos = v3_add(cam_pos, v3_scale(K_cam, -speed * (f32)is_key_down(KB_W)));
        cam_pos = v3_add(cam_pos, v3_scale(K_cam,  speed * (f32)is_key_down(KB_S)));

        cam_pos = v3_add(cam_pos, v3_scale(I_cam,  speed * (f32)is_key_down(KB_D)));
        cam_pos = v3_add(cam_pos, v3_scale(I_cam, -speed * (f32)is_key_down(KB_A)));
        
        cam_pos = v3_add(cam_pos, v3_scale(J_cam,  speed * (f32)is_key_down(KB_SPACE)));
        cam_pos = v3_add(cam_pos, v3_scale(J_cam, -speed * (f32)is_key_down(KB_LCTRL)));

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
         * Z should be between -1 and 1, so we need to divide by far plane - near plane. Use the vertex's Z coordinate instead of Q's Z coordinate
         * (Q is already projected and will have a constant Z, so we can't use that).
         * f : far plane dist
         *
         * Q_pz = (V_z - (-n)) * 2
         *        ----------------  -  1
         *           -f - (-n)
         *
         * Q_pz = (V_z + n) * 2
         *        -------------  -  1
         *           n - f
         *
         * Q_pz = V_z*2 + n*2
         *        -----------  -  1
         *           n - f
         *
         * Q_pz = V_z*2       n*2
         *        ------  +  -----  -  1
         *        n - f      n - f
         *
         * Q_pz =   2              n*2
         *        ----- * V_z  +  -----  -  1
         *        n - f           n - f
         *
         * https://www.desmos.com/calculator/frzetn7doc
         *       
         * Now, define as a matrix (keep in mind we will be dividing by the W component after matrix multiplication):
         *
         * | Q_px |   | n / C_w    0           0            0           |   |  V_x |   
         * | Q_py | = |   0      n / C_h       0            0           | * |  V_y |
         * | Q_pz |   |   0        0        2 / (n-f)    (n*2) / (n-f)  |   |  V_z |   
         * | Q_pw |   |   0        0          -1            0           |   | 1.0f |   
         *
         * | Q_px |   |         (n / C_w) * V_x         |
         * | Q_py | = |         (n / C_h) * V_y         |
         * | Q_pz |   | 2 / (n-f) * V_z + (n*2) / (n-f) |
         * | Q_pw |   |              -V_z               | <-- Will be dividing all the terms by -V_z
         * 
         * Dividing the depth by -V_z has the unfortunate consequence of reducing the NDC depth space.
         *
         */

        f32 aspect_ratio = get_screen_aspect_ratio();

        f32 n = 0.1f;
        f32 f = 1000.0f;
        f32 C_w = 0.125f;
        f32 C_h = C_w * aspect_ratio;
        mtx4x4 proj_mtx = {
            .m = {
                //    X         Y                Z                      W
                n / C_w,     0.0f,            0.0f,                  0.0f,
                   0.0f,  n / C_h,            0.0f,                  0.0f,
                   0.0f,     0.0f,  2.0f / (n - f),  (n * 2.0f) / (n - f),
                   0.0f,     0.0f,           -1.0f,                  0.0f,
            }
        };

        mtx4x4 camera_and_clip_mtx;
        mtx4x4_mul(&camera_and_clip_mtx, &proj_mtx, &world_to_cam_mtx);

        f32 pos[3] = {0};
        draw_Mesh_1M(terrain_mesh, &camera_and_clip_mtx, 1, pos + 0, pos + 1, pos + 2);

#if 0
        // Draw marching cube vertices
        {
            f32 circle_radius = 0.05f;
            v3 c_white = {.m={1.0f, 1.0f, 1.0f}};
            v3 c_black = {.m={0.0f, 0.0f, 0.0f}};

            for(u64_m i = 0; i < num_grid_points; i++)
            {
                debug_draw_sphere(&proj_mtx, &world_to_cam_mtx, grid_points[i], circle_radius, grid_points_filled[i] ? c_white : c_black);
            }
        }
#endif

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

}


