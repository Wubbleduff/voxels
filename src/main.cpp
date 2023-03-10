
#include <windows.h>
#include <stdio.h> // fopen
#include <immintrin.h>

#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "common.h"
#include "game_math.h"

/*
* Terrain scrolling in all directions
* Improve runtime of generating chunks
*/


// NOTE: This value must match the shader.
static const u32 SHADER_BUFFER_WIDTH = 1024*1024;
static const u32 NUM_VOXEL_DATA_FIELDS = 4;

static const u32 MAX_VOXELS = 50*1024*1024;
struct VoxelRenderData
{
    u32 num = 0;
    // | xxxx yyyy zzz |
    f32 pos[MAX_VOXELS*3];
    u32 color[MAX_VOXELS];

    static const u32 BATCH_SIZE = 10*1024;
};

struct Camera
{
    v3 pos;
    f32 rot_x;
    f32 rot_y;
    f32 vfov;
    f32 view_dist;
    f32 near_plane_dist;
    f32 ar; // w / h
};

struct GraphicsState
{
    GLFWwindow *window;
    char debug_info_log[1024];
    
    GLuint batch_voxel_shader_program;
    GLuint batch_voxel_vao;
    GLuint batch_voxel_vbo;
    GLuint batch_voxel_ebo;
    GLuint batch_voxel_ssbo;

    GLuint debug_line_shader_program;
    GLuint debug_line_vao;
    GLuint debug_line_vbo;

    Camera cam;
    Camera debug_cam;
    bool use_debug_cam;

    f32 voxel_vertices[32] =
    {
         0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f, // X
         0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f, // Y
        -0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f  // Z
    };
};
static GraphicsState *graphics_state = nullptr;

struct InputState
{
    s32 mouse_screen_x;
    s32 mouse_screen_y;
};
static InputState *input_state = nullptr;

static constexpr s32 CHUNK_DIM = 32;
static constexpr s32 CHUNK_POW = 5; // CHUNK_DIM = 2**CHUNK_POW
static constexpr s32 VIEW_DIST = 128;
static constexpr s32 LOADED_REGION_CHUNKS_DIM = (VIEW_DIST/CHUNK_DIM)*2;
struct ChunkBookmarkData
{
    static constexpr u32 MAX_CHUNKS = LOADED_REGION_CHUNKS_DIM*LOADED_REGION_CHUNKS_DIM*LOADED_REGION_CHUNKS_DIM;
    u32 num = 0;
    // | xxxx yyyy zzz |
    s32 pos[MAX_CHUNKS*3];
    u32 bookmark[MAX_CHUNKS+1];
};
struct VoxelData
{
    u32 num = 0;
    // | xxxx yyyy zzz |
    s32 pos[MAX_VOXELS*3];
    u32 color[MAX_VOXELS];
};

struct Timer
{
    Timer() { QueryPerformanceCounter(&start); }
    s64 time()
    {
        LARGE_INTEGER end;
        QueryPerformanceCounter(&end);
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        LARGE_INTEGER ms;
        ms.QuadPart = (end.QuadPart - start.QuadPart) * 1000;
        ms.QuadPart /= freq.QuadPart;
        return ms.QuadPart;
    }
    LARGE_INTEGER start;
};

struct ScopeTimer
{
    ScopeTimer(const char *name)
    {
        m_name = name;
        QueryPerformanceCounter(&start);
    }
    ~ScopeTimer()
    {
        LARGE_INTEGER end;
        QueryPerformanceCounter(&end);

        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);

        LARGE_INTEGER ms;
        ms.QuadPart = (end.QuadPart - start.QuadPart) * 1000;
        ms.QuadPart /= freq.QuadPart;
        ImGui::Text("TIME - %s: %ims", m_name, ms.QuadPart);
    }
    const char *m_name;
    LARGE_INTEGER start;
};
#define TIME_SCOPE(name) ScopeTimer _time_scope = ScopeTimer(name)



// TODO speed
mat4 view_m_world(const Camera* __restrict cam)
{
    mat4 y_rot = make_y_axis_rotation_matrix(-cam->rot_y);
    mat4 x_rot = make_x_axis_rotation_matrix(-cam->rot_x);
    mat4 result = x_rot * y_rot * make_translation_matrix(-cam->pos);
    return result;
}
mat4 clip_m_view(const Camera* __restrict cam)
{
    static f32 n = cam->near_plane_dist;
    static f32 f = cam->view_dist;
    f32 fov = deg_to_rad(cam->vfov);
    f32 r = -(f + n) / (f - n);
    f32 s = -(2.0f * n * f) / (f - n);
    mat4 result =
    {
         (f32)(1.0f / (tanf(fov / 2.0f) * cam->ar)) , 0.0f, 0.0f, 0.0f,
         0.0f, 1.0f / tanf(fov / 2.0f), 0.0f, 0.0f,
         0.0f, 0.0f, r, s,
         0.0f, 0.0f, -1.0f, 0.0f
    };
    return result;
}

// Outputs plan normals for right, top, left, bottom, near, far in that order.
static void get_frustum_planes(v3* __restrict out_normals, v3* __restrict out_points, const Camera* __restrict cam)
{
    mat4 cam_xform = view_m_world(&graphics_state->cam);
    v3 cam_p = graphics_state->cam.pos;
    v3 cam_i = v3(cam_xform[0][0], cam_xform[0][1], cam_xform[0][2]);
    v3 cam_j = v3(cam_xform[1][0], cam_xform[1][1], cam_xform[1][2]);
    v3 cam_k = v3(cam_xform[2][0], cam_xform[2][1], cam_xform[2][2]);

    f32 tan_fov = tan(deg_to_rad(graphics_state->cam.vfov)*0.5f);
    f32 near_hh = tan_fov * graphics_state->cam.near_plane_dist;
    f32 near_hw = near_hh * graphics_state->cam.ar;
    f32 far_hh = tan_fov * graphics_state->cam.view_dist;
    f32 far_hw = far_hh * graphics_state->cam.ar;

    v3 p_near = cam_p + -cam_k * graphics_state->cam.near_plane_dist;
    v3 p_far = cam_p + -cam_k * graphics_state->cam.view_dist;

    v3 frustum_v[] =
    {
        p_near - cam_i*near_hw - cam_j*near_hh, // bl 
        p_near + cam_i*near_hw - cam_j*near_hh, // br 
        p_near + cam_i*near_hw + cam_j*near_hh, // tr 
        p_near - cam_i*near_hw + cam_j*near_hh, // tl 

        p_far - cam_i*far_hw - cam_j*far_hh, // bl 
        p_far + cam_i*far_hw - cam_j*far_hh, // br 
        p_far + cam_i*far_hw + cam_j*far_hh, // tr 
        p_far - cam_i*far_hw + cam_j*far_hh, // tl 
    };
    // TODO Think harder about this...
    out_normals[0] = normalize(cross(frustum_v[6] - frustum_v[5], frustum_v[1] - frustum_v[5]));
    out_normals[1] = normalize(cross(frustum_v[7] - frustum_v[6], frustum_v[2] - frustum_v[6]));
    out_normals[2] = normalize(cross(frustum_v[4] - frustum_v[7], frustum_v[3] - frustum_v[7]));
    out_normals[3] = normalize(cross(frustum_v[5] - frustum_v[4], frustum_v[0] - frustum_v[4]));
    out_normals[4] =  cam_k;
    out_normals[5] = -cam_k;

    out_points[0] = cam_p;
    out_points[1] = cam_p;
    out_points[2] = cam_p;
    out_points[3] = cam_p;
    out_points[4] = p_near;
    out_points[5] = p_far;
}


static constexpr u32 left_pack_lut[] = 
{
    0x00000000, 0x00000000, 0x00000001, 0x00000010, 0x00000002, 0x00000020, 0x00000021, 0x00000210,
    0x00000003, 0x00000030, 0x00000031, 0x00000310, 0x00000032, 0x00000320, 0x00000321, 0x00003210,
    0x00000004, 0x00000040, 0x00000041, 0x00000410, 0x00000042, 0x00000420, 0x00000421, 0x00004210,
    0x00000043, 0x00000430, 0x00000431, 0x00004310, 0x00000432, 0x00004320, 0x00004321, 0x00043210,
    0x00000005, 0x00000050, 0x00000051, 0x00000510, 0x00000052, 0x00000520, 0x00000521, 0x00005210,
    0x00000053, 0x00000530, 0x00000531, 0x00005310, 0x00000532, 0x00005320, 0x00005321, 0x00053210,
    0x00000054, 0x00000540, 0x00000541, 0x00005410, 0x00000542, 0x00005420, 0x00005421, 0x00054210,
    0x00000543, 0x00005430, 0x00005431, 0x00054310, 0x00005432, 0x00054320, 0x00054321, 0x00543210,
    0x00000006, 0x00000060, 0x00000061, 0x00000610, 0x00000062, 0x00000620, 0x00000621, 0x00006210,
    0x00000063, 0x00000630, 0x00000631, 0x00006310, 0x00000632, 0x00006320, 0x00006321, 0x00063210,
    0x00000064, 0x00000640, 0x00000641, 0x00006410, 0x00000642, 0x00006420, 0x00006421, 0x00064210,
    0x00000643, 0x00006430, 0x00006431, 0x00064310, 0x00006432, 0x00064320, 0x00064321, 0x00643210,
    0x00000065, 0x00000650, 0x00000651, 0x00006510, 0x00000652, 0x00006520, 0x00006521, 0x00065210,
    0x00000653, 0x00006530, 0x00006531, 0x00065310, 0x00006532, 0x00065320, 0x00065321, 0x00653210,
    0x00000654, 0x00006540, 0x00006541, 0x00065410, 0x00006542, 0x00065420, 0x00065421, 0x00654210,
    0x00006543, 0x00065430, 0x00065431, 0x00654310, 0x00065432, 0x00654320, 0x00654321, 0x06543210,
    0x00000007, 0x00000070, 0x00000071, 0x00000710, 0x00000072, 0x00000720, 0x00000721, 0x00007210,
    0x00000073, 0x00000730, 0x00000731, 0x00007310, 0x00000732, 0x00007320, 0x00007321, 0x00073210,
    0x00000074, 0x00000740, 0x00000741, 0x00007410, 0x00000742, 0x00007420, 0x00007421, 0x00074210,
    0x00000743, 0x00007430, 0x00007431, 0x00074310, 0x00007432, 0x00074320, 0x00074321, 0x00743210,
    0x00000075, 0x00000750, 0x00000751, 0x00007510, 0x00000752, 0x00007520, 0x00007521, 0x00075210,
    0x00000753, 0x00007530, 0x00007531, 0x00075310, 0x00007532, 0x00075320, 0x00075321, 0x00753210,
    0x00000754, 0x00007540, 0x00007541, 0x00075410, 0x00007542, 0x00075420, 0x00075421, 0x00754210,
    0x00007543, 0x00075430, 0x00075431, 0x00754310, 0x00075432, 0x00754320, 0x00754321, 0x07543210,
    0x00000076, 0x00000760, 0x00000761, 0x00007610, 0x00000762, 0x00007620, 0x00007621, 0x00076210,
    0x00000763, 0x00007630, 0x00007631, 0x00076310, 0x00007632, 0x00076320, 0x00076321, 0x00763210,
    0x00000764, 0x00007640, 0x00007641, 0x00076410, 0x00007642, 0x00076420, 0x00076421, 0x00764210,
    0x00007643, 0x00076430, 0x00076431, 0x00764310, 0x00076432, 0x00764320, 0x00764321, 0x07643210,
    0x00000765, 0x00007650, 0x00007651, 0x00076510, 0x00007652, 0x00076520, 0x00076521, 0x00765210,
    0x00007653, 0x00076530, 0x00076531, 0x00765310, 0x00076532, 0x00765320, 0x00765321, 0x07653210,
    0x00007654, 0x00076540, 0x00076541, 0x00765410, 0x00076542, 0x00765420, 0x00765421, 0x07654210,
    0x00076543, 0x00765430, 0x00765431, 0x07654310, 0x00765432, 0x07654320, 0x07654321, 0x76543210 
};
static __m256 left_pack(__m256 a, u32 mask)
{
    __m256i shufmask = _mm256_srlv_epi32(_mm256_set1_epi32(left_pack_lut[mask]), _mm256_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28));
    return _mm256_permutevar8x32_ps(a, shufmask);
}
static __m256i left_pack(__m256i a, u32 mask)
{
    __m256i shufmask = _mm256_srlv_epi32(_mm256_set1_epi32(left_pack_lut[mask]), _mm256_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28));
    return _mm256_permutevar8x32_epi32(a, shufmask);
}

void camera_cull(
    u32* __restrict out_num_voxels,
    f32* __restrict out_voxel_pos,
    u32* __restrict out_voxel_color,

    const u32 in_num_voxels,
    const s32* __restrict in_voxel_pos,
    const u32* __restrict in_voxel_color,
    const Camera* __restrict cam)
{
#if 0
    // Pass through
    voxel_render_data->num = 0;
    for(u32 i = 0; i < all_voxel_data->num; i++)
    {
        voxel_render_data->x[voxel_render_data->num] = all_voxel_data->x[i];
        voxel_render_data->y[voxel_render_data->num] = all_voxel_data->y[i];
        voxel_render_data->z[voxel_render_data->num] = all_voxel_data->z[i];
        voxel_render_data->color[voxel_render_data->num] = all_voxel_data->color[i];
        voxel_render_data->num++;
    }
#elif 0
    // Reference
    const mat4 clip_mat = clip_m_view(&graphics_state->cam) * view_m_world(&graphics_state->cam);
    f32 * __restrict vv = graphics_state->voxel_vertices;
    voxel_render_data->num = 0;
    for(u32 i = 0; i < all_voxel_data->num; i++)
    {
        f32 voxel_x = all_voxel_data->x[i];
        f32 voxel_y = all_voxel_data->y[i];
        f32 voxel_z = all_voxel_data->z[i];
        u32 voxel_color = all_voxel_data->color[i];

        bool cull = true;
        for(u32 vi = 0; vi < 8; vi++)
        {
            v4 v = v4(voxel_x + vv[vi], voxel_y + vv[8+vi], voxel_z + vv[16+vi], 1.0f);
            v4 v_clip = clip_mat * v;
            if((v_clip.x > -v_clip.w) && (v_clip.x < v_clip.w) && 
               (v_clip.y > -v_clip.w) && (v_clip.y < v_clip.w) && 
               (v_clip.z > -v_clip.w) && (v_clip.z < v_clip.w))
            {
                cull = false;
            }
        }

        if(!cull)
        {
            voxel_render_data->x[voxel_render_data->num] = voxel_x;
            voxel_render_data->y[voxel_render_data->num] = voxel_y;
            voxel_render_data->z[voxel_render_data->num] = voxel_z;
            voxel_render_data->color[voxel_render_data->num] = voxel_color;
            voxel_render_data->num++;
        }
    }
#else
    // Check distance from voxel centroid to each plane and compare against 1.0f.
    // This is effectively the same as a sphere-plane distance check.
    // If P is a point on the plane, Q is the voxel centroid, and n is the plane normal, dist check is:
    // (Q - P) * n < sqrt(0.5^2 + 0.5^2) ->
    // Q*n - P*n - sqrt(0.5^2 + 0.5^2)
    // We can precomupte p_dot_n = P*n + sqrt(0.5^2 + 0.5^2) outside the per-voxel loop.
    // Inside the loop, we just need to compute Q*n - p_dot_n. The resulting dist will
    // have the sign bit set if < 0. Then we can use movemask without needing an extra
    // cmp.

    v3 plane_normals[6];
    v3 plane_points[6];
    get_frustum_planes(plane_normals, plane_points, cam);

    __m256 p_dot_n[] =
    {
        // Add 0.1f for padding. Otherwise we see overly optimistic culling.
        _mm256_set1_ps(dot(plane_points[0], plane_normals[0]) + 0.707107f+0.1f),
        _mm256_set1_ps(dot(plane_points[1], plane_normals[1]) + 0.707107f+0.1f),
        _mm256_set1_ps(dot(plane_points[2], plane_normals[2]) + 0.707107f+0.1f),
        _mm256_set1_ps(dot(plane_points[3], plane_normals[3]) + 0.707107f+0.1f),
        _mm256_set1_ps(dot(plane_points[4], plane_normals[4]) + 0.707107f+0.1f),
        _mm256_set1_ps(dot(plane_points[5], plane_normals[5]) + 0.707107f+0.1f)
    };

    u32 num_voxels = 0;
    for(u32 i = 0; i < in_num_voxels; i += 8)
    {
        __m256 vx = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(in_voxel_pos + MAX_VOXELS*0 + i)));
        __m256 vy = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(in_voxel_pos + MAX_VOXELS*1 + i)));
        __m256 vz = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(in_voxel_pos + MAX_VOXELS*2 + i)));
        __m256i vc = _mm256_loadu_si256((__m256i*)(in_voxel_color + i));

        __m256 mask = _mm256_cvtepi32_ps(_mm256_set1_epi8(0xFF));
        for(u32 ni = 0; ni < 6; ni++)
        {
            __m256 d = _mm256_fmadd_ps(vx, _mm256_broadcast_ss(&plane_normals[ni].x),
                       _mm256_fmadd_ps(vy, _mm256_broadcast_ss(&plane_normals[ni].y),
                       _mm256_mul_ps(vz, _mm256_broadcast_ss(&plane_normals[ni].z))));

            __m256 dist = _mm256_sub_ps(d, p_dot_n[ni]);
            mask = _mm256_and_ps(mask, dist);
        }

        u32 pack_index = _mm256_movemask_ps(mask);
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*0 + num_voxels,   left_pack(vx, pack_index));
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*1 + num_voxels,   left_pack(vy, pack_index));
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*2 + num_voxels,   left_pack(vz, pack_index));
        _mm256_storeu_si256((__m256i*)(out_voxel_color + num_voxels), left_pack(vc, pack_index));
        num_voxels += _mm_popcnt_u32(pack_index);
    }
    *out_num_voxels = num_voxels;
#endif
}


void glfw_error_fn(s32 error, const char *message)
{
    //printf("%i: %s\n", error, message);
    assert(false);
}

static void check_gl_errors(const char *desc)
{
    GLint error = glGetError();
    if(error)
    {
        glfw_error_fn(error, desc);
    }
}

static char *read_file_into_string(const char *path)
{
    FILE *file = fopen(path, "rb");
    if(file == nullptr) return nullptr;
    fseek(file, 0L, SEEK_END);
    u32 size = ftell(file);
    fseek(file, 0L, SEEK_SET);
    char *buffer = new char[size + 1]();
    fread(buffer, size, 1, file);
    buffer[size] = '\0';
    return buffer;
}

static bool is_newline(char *c)
{
    return (c[0] == '\n' || c[0] == '\r' || (c[0] == '\r' && c[1] == '\n'));
}

static bool read_shader_file(const char *path, char **vert_source, char **geom_source, char **frag_source)
{
    *vert_source = nullptr;
    *geom_source = nullptr;
    *frag_source = nullptr;
    
    char *file = read_file_into_string(path);
    if(file == nullptr) return false;
    
    u32 current_tag_length = 0;
    char *current_tag = nullptr;
    char *current_source_start = nullptr;
    
    char *character = file;
    while(*character != '\0')
    {
        if(*character == '@')
        {
            // Finish reading current shader source
            if(current_tag == nullptr)
            {
            }
            else if(strncmp("vertex", current_tag, current_tag_length) == 0)
            {
                *vert_source = current_source_start;
            }
            else if(strncmp("geometry", current_tag, current_tag_length) == 0)
            {
                *geom_source = current_source_start;
            }
            else if(strncmp("fragment", current_tag, current_tag_length) == 0)
            {
                *frag_source = current_source_start;
            }
            
            
            // Null terminate previous shader string
            *character = '\0';
            
            // Read tag
            character++;
            char *tag = character;
            
            // Move past tag
            while(!is_newline(character))
            {
                character++;
            }
            char *one_past_end_tag = character;
            while(is_newline(character)) character++;
            
            current_tag_length = one_past_end_tag - tag;
            current_tag = tag;
            current_source_start = character;
        }
        else
        {
            character++;
        }
    }
    
    // Finish reading current shader source
    if(current_tag == nullptr)
    {
    }
    else if(strncmp("vertex", current_tag, current_tag_length) == 0)
    {
        *vert_source = current_source_start;
    }
    else if(strncmp("geometry", current_tag, current_tag_length) == 0)
    {
        *geom_source = current_source_start;
    }
    else if(strncmp("fragment", current_tag, current_tag_length) == 0)
    {
        *frag_source = current_source_start;
    }
    
    return true;
}

GLuint make_shader_from_string(const char *vert_source, const char *frag_source)
{
    s32 success;
    const u32 info_log_size = sizeof(graphics_state->debug_info_log);
    
    u32 vert_shader;
    vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vert_source, NULL);
    glCompileShader(vert_shader);
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(vert_shader, info_log_size, NULL, graphics_state->debug_info_log);
        return -1;
    }
    
    u32 frag_shader;
    frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(frag_shader, 1, &frag_source, NULL);
    glCompileShader(frag_shader);
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(frag_shader, info_log_size, NULL, graphics_state->debug_info_log);
        return -1;
    }
    
    check_gl_errors("compiling shaders");
    
    GLuint program = glCreateProgram();
    check_gl_errors("making program");
    
    glAttachShader(program, vert_shader);
    glAttachShader(program, frag_shader);
    glLinkProgram(program);
    check_gl_errors("linking program");
    
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success)
    {
        glGetProgramInfoLog(program, info_log_size, NULL, graphics_state->debug_info_log);
        return -1;
    }
    
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader); 
    check_gl_errors("deleting shaders");
    
    return program;
}

GLuint make_shader_from_file(const char *shader_path)
{
    s32 success;
    const u32 info_log_size = sizeof(graphics_state->debug_info_log);
    
    char *vert_source = nullptr;
    char *geom_source = nullptr;
    char *frag_source = nullptr;
    bool read = read_shader_file(shader_path, &vert_source, &geom_source, &frag_source);
    if(!read)
    {
        //fprintf(stderr, "Could not open shader file %s\n", shader_path);
        return -1;
    }

    return make_shader_from_string(vert_source, frag_source);
}

static void reshape(GLFWwindow *window, s32 width, s32 height)
{
    glViewport(0, 0, (GLint)width, (GLint)height);
    graphics_state->cam.ar = (f32)width / height;
    graphics_state->debug_cam.ar = (f32)width / height;
}

#if 0
// Right now I'm just defining "chunk" as a MxNxO group of voxels
static void terrain_chunk(VoxelData *out,
        s32 offset_x, s32 offset_z,
        s32 width, s32 depth,
        f32 noise_scale_x, f32 noise_scale_y, f32 noise_scale_z)
{
    auto noise = [noise_scale_x, noise_scale_y, noise_scale_z](s32 x_, s32 z_)
    {
        f32 n = (perlin_noise(x_*noise_scale_x, z_*noise_scale_z)*0.5f + 0.5f)*noise_scale_y;
        n += (perlin_noise(x_*noise_scale_x*0.1f, z_*noise_scale_z*0.1f)*0.5f + 0.5f)*noise_scale_y*10.0f;
        return n;
    };
    auto pick_color = [](s32 y)
    {
        if(y < 65) return (u32)0x424242FF;
        else if(y < 480) return (u32)0x519916FF;
        else return (u32)0xCCCCCCFF;
    };

    u32 num = 0;
    for(s32 z = offset_z; z < depth + offset_z; z++)
    {
        for(s32 x = offset_x; x < width + offset_x; x++)
        {
            s32 y = (s32)noise(x, z);

            // Check for holes in the ground
            {
                s32 y0 = (s32)noise(x - 1, z);
                s32 y1 = (s32)noise(x + 1, z);
                s32 y2 = (s32)noise(x, z - 1);
                s32 y3 = (s32)noise(x, z + 1);
                s32 max_d = max(y - y0, max(y - y1, max(y - y2, y - y3)));
                for(s32 i = 0; i < max_d - 1; i++)
                {
                    s32 yi = y - 1 - i;
                    out->pos[MAX_VOXELS*0 + num] = x;
                    out->pos[MAX_VOXELS*1 + num] = yi;
                    out->pos[MAX_VOXELS*2 + num] = z;
                    out->color[num] = pick_color(yi);
                    num++;
                }
            }

            out->pos[MAX_VOXELS*0 + num] = x;
            out->pos[MAX_VOXELS*1 + num] = y;
            out->pos[MAX_VOXELS*2 + num] = z;
            out->color[num] = pick_color(y);
            num++;
        }
    }
}
#endif

static void draw_line(v3 a, v3 b, v3 c)
{
    glUseProgram(graphics_state->debug_line_shader_program);
    check_gl_errors("use program");

    Camera *current_cam;
    if(graphics_state->use_debug_cam) current_cam = &graphics_state->debug_cam;
    else current_cam = &graphics_state->cam;

    mat4 m = clip_m_view(current_cam) * view_m_world(current_cam);
    GLint loc = glGetUniformLocation(graphics_state->debug_line_shader_program, "mvp");
    glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
    if(loc == -1) assert(false);

    loc = glGetUniformLocation(graphics_state->debug_line_shader_program, "color");
    glUniform3f(loc, c.r, c.g, c.b);
    if(loc == -1) assert(false);

    glBindVertexArray(graphics_state->debug_line_vao);
    glBindBuffer(GL_ARRAY_BUFFER, graphics_state->debug_line_vbo);

    f32 v[] =
    {
        a.x, a.y, a.z,
        b.x, b.y, b.z
    };
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(v), v);

    glDrawArrays(GL_LINES, 0, 2);
    check_gl_errors("draw");
}


// Chunks
// X
// Y
// Z
// # voxels
// voxels bookmark

f32 terrain_noise_3d(s32 x, s32 y, s32 z)
{
    return pnoise(x * 0.01f, y * 0.1f, z * 0.01f);
}

void gen_chunk(
    ChunkBookmarkData* __restrict out_chunk_data,
    VoxelData* __restrict out_voxel_data,
    const s32 chunk_x,
    const s32 chunk_y,
    const s32 chunk_z
)
{
    constexpr f32 y_bias = -0.1f;

//#define DO_CHUNK_SAMPLE_OPT
#ifdef DO_CHUNK_SAMPLE_OPT
    {
        s32 x0 = chunk_x * CHUNK_DIM;
        s32 x1 = chunk_x * CHUNK_DIM + CHUNK_DIM;
        s32 y0 = chunk_y * CHUNK_DIM;
        s32 y1 = chunk_y * CHUNK_DIM + CHUNK_DIM;
        s32 z0 = chunk_z * CHUNK_DIM;
        s32 z1 = chunk_z * CHUNK_DIM + CHUNK_DIM;
        s32 samples[] =
        {
            x0, y0, z0,
            x0, y0, z1,
            x0, y1, z0,
            x0, y1, z1,
            x1, y0, z0,
            x1, y0, z1,
            x1, y1, z0,
            x1, y1, z1,
        };

        bool all_in  = true;
        bool all_out = true;
        for(u32 i = 0; i < 8; i++)
        {
            s32 x = samples[i*3 + 0];
            s32 y = samples[i*3 + 1];
            s32 z = samples[i*3 + 2];
            f32 n = terrain_noise_3d(x, y, z) - y*y_bias;
            if(n <  0.0f) all_out = false;
            if(n >= 0.0f) all_in  = false;
        }
        if(all_in || all_out)
        {
            return;
        }
    }
#endif

    u32* __restrict bookmark = out_chunk_data->bookmark;
    s32* __restrict voxels_x = &(out_voxel_data->pos[MAX_VOXELS*0]);
    s32* __restrict voxels_y = &(out_voxel_data->pos[MAX_VOXELS*1]);
    s32* __restrict voxels_z = &(out_voxel_data->pos[MAX_VOXELS*2]);
    u32* __restrict voxels_color = out_voxel_data->color;

    const u32 n_c = out_chunk_data->num;
    out_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*0 + n_c] = chunk_x;
    out_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*1 + n_c] = chunk_y;
    out_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*2 + n_c] = chunk_z;

    if(n_c == 0)
    {
        bookmark[0] = 0;
    }

    u32 base_voxels = bookmark[n_c];
    u32 num_voxels = 0;
    for(s32 x = chunk_x * CHUNK_DIM; x < chunk_x * CHUNK_DIM + CHUNK_DIM; x++)
    {
        for(s32 y = chunk_y * CHUNK_DIM; y < chunk_y * CHUNK_DIM + CHUNK_DIM; y++)
        {
            for(s32 z = chunk_z * CHUNK_DIM; z < chunk_z * CHUNK_DIM + CHUNK_DIM; z++)
            {
                f32 n = terrain_noise_3d(x, y, z) - y*y_bias;
                if(n < 0.0f)
                {
                    voxels_x[base_voxels + num_voxels] = x;
                    voxels_y[base_voxels + num_voxels] = y;
                    voxels_z[base_voxels + num_voxels] = z;
                    voxels_color[base_voxels + num_voxels] = 0x00AA00FF;
                    ++num_voxels;
                }
            }
        }
    }
    bookmark[n_c + 1] = base_voxels + num_voxels;
    out_chunk_data->num = n_c + 1;
}

// TODO Fix issue where we can jump multiple chunks at a time.... Right now we only assume 1 chunk can be jumped at a time.
bool maybe_gen_new_terrain(
    ChunkBookmarkData* __restrict dst_chunk_data,
    VoxelData* __restrict dst_voxel_data,

    const ChunkBookmarkData* __restrict src_chunk_data,
    const VoxelData* __restrict src_voxel_data,

    const s32 old_chunk_x,
    const s32 old_chunk_y,
    const s32 old_chunk_z,
    const s32 new_chunk_x,
    const s32 new_chunk_y,
    const s32 new_chunk_z
)
{
    bool did_generation = false;

    s32* __restrict dst_chunks_x = &dst_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*0];
    s32* __restrict dst_chunks_y = &dst_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*1];
    s32* __restrict dst_chunks_z = &dst_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*2];
    s32* __restrict dst_voxels_x = &dst_voxel_data->pos[MAX_VOXELS*0];
    s32* __restrict dst_voxels_y = &dst_voxel_data->pos[MAX_VOXELS*1];
    s32* __restrict dst_voxels_z = &dst_voxel_data->pos[MAX_VOXELS*2];
    u32* __restrict dst_voxels_color = dst_voxel_data->color;

    const s32* __restrict src_chunks_x = &src_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*0];
    const s32* __restrict src_chunks_y = &src_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*1];
    const s32* __restrict src_chunks_z = &src_chunk_data->pos[ChunkBookmarkData::MAX_CHUNKS*2];
    const s32* __restrict src_voxels_x = &src_voxel_data->pos[MAX_VOXELS*0];
    const s32* __restrict src_voxels_y = &src_voxel_data->pos[MAX_VOXELS*1];
    const s32* __restrict src_voxels_z = &src_voxel_data->pos[MAX_VOXELS*2];
    const u32* __restrict src_voxels_color = src_voxel_data->color;

    if(old_chunk_x != new_chunk_x)
    {
        // NOTE: LOADED_REGION_CHUNKS_DIM is assumed even.
        s32 delete_chunk_coord;
        s32 new_chunk_coord;
        if(old_chunk_x < new_chunk_x)
        {
            delete_chunk_coord = old_chunk_x - LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = old_chunk_x + LOADED_REGION_CHUNKS_DIM/2;
        }
        else
        {
            delete_chunk_coord = new_chunk_x + LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = new_chunk_x - LOADED_REGION_CHUNKS_DIM/2;
        }

        u32 dst_chunks_idx = 0;
        u32 dst_voxels_idx = 0;
        dst_chunk_data->bookmark[0] = 0;
        for(u32 src_chunk_idx = 0; src_chunk_idx < src_chunk_data->num; src_chunk_idx++)
        {
            if(src_chunks_x[src_chunk_idx] != delete_chunk_coord)
            {
                dst_chunks_x[dst_chunks_idx] = src_chunks_x[src_chunk_idx];
                dst_chunks_y[dst_chunks_idx] = src_chunks_y[src_chunk_idx];
                dst_chunks_z[dst_chunks_idx] = src_chunks_z[src_chunk_idx];

                // Left pack voxels
                const u32 src_num_voxels = src_chunk_data->bookmark[src_chunk_idx + 1] - src_chunk_data->bookmark[src_chunk_idx];
                const u32 src_voxels_idx = src_chunk_data->bookmark[src_chunk_idx];
                memcpy(dst_voxels_x + dst_voxels_idx, src_voxels_x + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_y + dst_voxels_idx, src_voxels_y + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_z + dst_voxels_idx, src_voxels_z + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_color + dst_voxels_idx, src_voxels_color + src_voxels_idx, src_num_voxels * sizeof(u32));

                dst_chunk_data->bookmark[dst_chunks_idx+1] = dst_voxels_idx + src_num_voxels;
                dst_voxels_idx += src_num_voxels;
                dst_chunks_idx++;
            }
        }
        assert(dst_chunks_idx <= ChunkBookmarkData::MAX_CHUNKS); // TODO
        dst_chunk_data->num = dst_chunks_idx;
        dst_chunk_data->bookmark[dst_chunk_data->num] = dst_voxels_idx;


        // Gen new terrain
        s32 region_chunks_bl_y = new_chunk_y - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_bl_z = new_chunk_z - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_y = new_chunk_y + LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_z = new_chunk_z + LOADED_REGION_CHUNKS_DIM/2;
        for(s32 chunk_y = region_chunks_bl_y; chunk_y < region_chunks_tr_y; chunk_y++)
        {
            for(s32 chunk_z = region_chunks_bl_z; chunk_z < region_chunks_tr_z; chunk_z++)
            {
                gen_chunk(
                    dst_chunk_data,
                    dst_voxel_data,
                    new_chunk_coord, chunk_y, chunk_z
                );
            }
        }
        dst_voxel_data->num = dst_chunk_data->bookmark[dst_chunk_data->num];
        did_generation = true;
    }

    if(old_chunk_y != new_chunk_y)
    {
        // NOTE: LOADED_REGION_CHUNKS_DIM is assumed even.
        s32 delete_chunk_coord;
        s32 new_chunk_coord;
        if(old_chunk_y < new_chunk_y)
        {
            delete_chunk_coord = old_chunk_y - LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = old_chunk_y + LOADED_REGION_CHUNKS_DIM/2;
        }
        else
        {
            delete_chunk_coord = new_chunk_y + LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = new_chunk_y - LOADED_REGION_CHUNKS_DIM/2;
        }

        u32 dst_chunks_idx = 0;
        u32 dst_voxels_idx = 0;
        dst_chunk_data->bookmark[0] = 0;
        for(u32 src_chunk_idx = 0; src_chunk_idx < src_chunk_data->num; src_chunk_idx++)
        {
            if(src_chunks_y[src_chunk_idx] != delete_chunk_coord)
            {
                dst_chunks_x[dst_chunks_idx] = src_chunks_x[src_chunk_idx];
                dst_chunks_y[dst_chunks_idx] = src_chunks_y[src_chunk_idx];
                dst_chunks_z[dst_chunks_idx] = src_chunks_z[src_chunk_idx];

                // Left pack voxels
                const u32 src_num_voxels = src_chunk_data->bookmark[src_chunk_idx + 1] - src_chunk_data->bookmark[src_chunk_idx];
                const u32 src_voxels_idx = src_chunk_data->bookmark[src_chunk_idx];
                memcpy(dst_voxels_x + dst_voxels_idx, src_voxels_x + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_y + dst_voxels_idx, src_voxels_y + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_z + dst_voxels_idx, src_voxels_z + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_color + dst_voxels_idx, src_voxels_color + src_voxels_idx, src_num_voxels * sizeof(u32));

                dst_chunk_data->bookmark[dst_chunks_idx+1] = dst_voxels_idx + src_num_voxels;
                dst_voxels_idx += src_num_voxels;
                dst_chunks_idx++;
            }
        }
        assert(dst_chunks_idx <= ChunkBookmarkData::MAX_CHUNKS); // TODO
        dst_chunk_data->num = dst_chunks_idx;
        dst_chunk_data->bookmark[dst_chunk_data->num] = dst_voxels_idx;

        // Gen new terrain
        s32 region_chunks_bl_x = new_chunk_x - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_bl_z = new_chunk_z - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_x = new_chunk_x + LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_z = new_chunk_z + LOADED_REGION_CHUNKS_DIM/2;
        for(s32 chunk_x = region_chunks_bl_x; chunk_x < region_chunks_tr_x; chunk_x++)
        {
            for(s32 chunk_z = region_chunks_bl_z; chunk_z < region_chunks_tr_z; chunk_z++)
            {
                gen_chunk(
                    dst_chunk_data,
                    dst_voxel_data,
                    chunk_x, new_chunk_coord, chunk_z
                );
            }
        }
        dst_voxel_data->num = dst_chunk_data->bookmark[dst_chunk_data->num];
        did_generation = true;
    }

    if(old_chunk_z != new_chunk_z)
    {
        // NOTE: LOADED_REGION_CHUNKS_DIM is assumed even.
        s32 delete_chunk_coord;
        s32 new_chunk_coord;
        if(old_chunk_z < new_chunk_z)
        {
            delete_chunk_coord = old_chunk_z - LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = old_chunk_z + LOADED_REGION_CHUNKS_DIM/2;
        }
        else
        {
            delete_chunk_coord = new_chunk_z + LOADED_REGION_CHUNKS_DIM/2;
            new_chunk_coord    = new_chunk_z - LOADED_REGION_CHUNKS_DIM/2;
        }

        u32 dst_chunks_idx = 0;
        u32 dst_voxels_idx = 0;
        dst_chunk_data->bookmark[0] = 0;
        for(u32 src_chunk_idx = 0; src_chunk_idx < src_chunk_data->num; src_chunk_idx++)
        {
            if(src_chunks_z[src_chunk_idx] != delete_chunk_coord)
            {
                dst_chunks_x[dst_chunks_idx] = src_chunks_x[src_chunk_idx];
                dst_chunks_y[dst_chunks_idx] = src_chunks_y[src_chunk_idx];
                dst_chunks_z[dst_chunks_idx] = src_chunks_z[src_chunk_idx];

                // Left pack voxels
                const u32 src_num_voxels = src_chunk_data->bookmark[src_chunk_idx + 1] - src_chunk_data->bookmark[src_chunk_idx];
                const u32 src_voxels_idx = src_chunk_data->bookmark[src_chunk_idx];
                memcpy(dst_voxels_x + dst_voxels_idx, src_voxels_x + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_y + dst_voxels_idx, src_voxels_y + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_z + dst_voxels_idx, src_voxels_z + src_voxels_idx, src_num_voxels * sizeof(f32));
                memcpy(dst_voxels_color + dst_voxels_idx, src_voxels_color + src_voxels_idx, src_num_voxels * sizeof(u32));

                dst_chunk_data->bookmark[dst_chunks_idx+1] = dst_voxels_idx + src_num_voxels;
                dst_voxels_idx += src_num_voxels;
                dst_chunks_idx++;
            }
        }
        assert(dst_chunks_idx <= ChunkBookmarkData::MAX_CHUNKS); // TODO
        dst_chunk_data->num = dst_chunks_idx;
        dst_chunk_data->bookmark[dst_chunk_data->num] = dst_voxels_idx;

        // Gen new terrain
        s32 region_chunks_bl_x = new_chunk_x - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_bl_y = new_chunk_y - LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_x = new_chunk_x + LOADED_REGION_CHUNKS_DIM/2;
        s32 region_chunks_tr_y = new_chunk_y + LOADED_REGION_CHUNKS_DIM/2;
        for(s32 chunk_x = region_chunks_bl_x; chunk_x < region_chunks_tr_x; chunk_x++)
        {
            for(s32 chunk_y = region_chunks_bl_y; chunk_y < region_chunks_tr_y; chunk_y++)
            {
                gen_chunk(
                    dst_chunk_data,
                    dst_voxel_data,
                    chunk_x, chunk_y, new_chunk_coord
                );
            }
        }
        dst_voxel_data->num = dst_chunk_data->bookmark[dst_chunk_data->num];
        did_generation = true;
    }

    return did_generation;
}



INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR lpCmdLine, INT nCmdShow)
{
    // GLFW
    graphics_state = new GraphicsState();
    {
        glfwSetErrorCallback(glfw_error_fn);
        if(!glfwInit())
        {
            //fprintf(stderr, "Failed to initialize GLFW\n");
            return 1;
        }
        
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_SAMPLES, 4); // TODO Not working?

        u32 window_width = 1920 * 1.0f;
        u32 window_height = 1080 * 1.0f;
        graphics_state->window = glfwCreateWindow(window_width, window_height, "My Game", NULL, NULL);
        if(!graphics_state->window)
        {
            //fprintf(stderr, "Failed to open GLFW window\n");
            glfwTerminate();
            return 1;
        }
        
        glfwSetFramebufferSizeCallback(graphics_state->window, reshape);
        //glfwSetKeyCallback(graphics_state->window, key);
        
        glfwMakeContextCurrent(graphics_state->window);
        
        if(!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
        {
            //fprintf(stderr, "Failed to initialize OpenGL context\n");
            return 1;
        }
        
        glfwSwapInterval(1);

        glEnable(GL_MULTISAMPLE); // TODO Not working?
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CCW);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        
        s32 framebuffer_width;
        s32 framebuffer_height;
        glfwGetFramebufferSize(graphics_state->window, &framebuffer_width, &framebuffer_height);
        reshape(graphics_state->window, framebuffer_width, framebuffer_height);

        graphics_state->cam.pos = v3(0.0f, 0.0f, -5.0f);
        graphics_state->cam.rot_x = 0.0f;
        graphics_state->cam.rot_y = M_PI;
        graphics_state->cam.vfov = 90.0f;
        graphics_state->cam.view_dist = 5000.0f;
        graphics_state->cam.near_plane_dist = 1.0f;
        graphics_state->debug_cam = graphics_state->cam;

        {
            glGenVertexArrays(1, &graphics_state->batch_voxel_vao);
            glBindVertexArray(graphics_state->batch_voxel_vao);
            check_gl_errors("vao");
            glGenBuffers(1, &graphics_state->batch_voxel_vbo);
            check_gl_errors("vbo");
            static u32 voxel_indices[] = 
            {
                4, 0, 3,
                4, 3, 7,
                2, 6, 7,
                2, 7, 3,
                5, 2, 1,
                5, 6, 2,
                4, 1, 0,
                4, 5, 1,
                7, 5, 4,
                7, 6, 5,
                0, 1, 2,
                0, 2, 3
            };

            glBindBuffer(GL_ARRAY_BUFFER, graphics_state->batch_voxel_vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(graphics_state->voxel_vertices), graphics_state->voxel_vertices, GL_STATIC_DRAW);
            check_gl_errors("vbo voxel vertices");

            glGenBuffers(1, &graphics_state->batch_voxel_ebo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_state->batch_voxel_ebo);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(voxel_indices), voxel_indices, GL_STATIC_DRAW);

            glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void *)0);
            glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (void *)(8*sizeof(f32)));
            glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, (void *)(16*sizeof(f32)));
            glEnableVertexAttribArray(0);
            glEnableVertexAttribArray(1);
            glEnableVertexAttribArray(2);
            check_gl_errors("vertex attrib pointer");

            // https://www.khronos.org/opengl/wiki/Shader_Storage_Buffer_Object
            const u32 voxel_buffer_size = SHADER_BUFFER_WIDTH * NUM_VOXEL_DATA_FIELDS * 4;
            glGenBuffers(1, &graphics_state->batch_voxel_ssbo);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, graphics_state->batch_voxel_ssbo);
            glBufferData(GL_SHADER_STORAGE_BUFFER, voxel_buffer_size, nullptr, GL_DYNAMIC_DRAW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, graphics_state->batch_voxel_ssbo);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            check_gl_errors("ssbo");

            graphics_state->batch_voxel_shader_program = make_shader_from_file("assets/shaders/batch_voxel.shader");
            {
                FILE *file = fopen("log.txt", "wt");
                fprintf(file, graphics_state->debug_info_log);
                fclose(file);
            }
        }

        {
            glGenVertexArrays(1, &graphics_state->debug_line_vao);
            glBindVertexArray(graphics_state->debug_line_vao);
            check_gl_errors("vao");
            glGenBuffers(1, &graphics_state->debug_line_vbo);
            check_gl_errors("vbo");
            glBindBuffer(GL_ARRAY_BUFFER, graphics_state->debug_line_vbo);
            f32 v[6] = {};
            glBufferData(GL_ARRAY_BUFFER, sizeof(v), v, GL_STATIC_DRAW);
            check_gl_errors("vbo voxel vertices");

            glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 3*sizeof(f32), (void *)0);               // x
            glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 3*sizeof(f32), (void *)(1*sizeof(f32))); // y
            glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 3*sizeof(f32), (void *)(2*sizeof(f32))); // z
            glEnableVertexAttribArray(0);
            glEnableVertexAttribArray(1);
            glEnableVertexAttribArray(2);
            check_gl_errors("vertex attrib pointer");

            graphics_state->debug_line_shader_program = make_shader_from_string(
                "#version 440 core\n"
                "layout (location = 0) in float a_x;"
                "layout (location = 1) in float a_y;"
                "layout (location = 2) in float a_z;"
                "uniform mat4 mvp;"
                "void main()"
                "{"
                "    gl_Position = mvp * vec4(a_x, a_y, a_z, 1.0f);"
                "};",

                "#version 440 core\n"
                "uniform vec3 color;"
                "out vec4 frag_color;"
                "void main()"
                "{"
                "    frag_color = vec4(color, 1);"
                "}"
            );
            {
                FILE *file = fopen("log.txt", "wt");
                fprintf(file, graphics_state->debug_info_log);
                fclose(file);
            }
        }
    }

    // ImGui
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplOpenGL3_Init("#version 440 core");
        ImGui_ImplGlfw_InitForOpenGL(graphics_state->window, true);
    }

    input_state = new InputState();

    ChunkBookmarkData * __restrict chunk_data0 = new ChunkBookmarkData;
    ChunkBookmarkData * __restrict chunk_data1 = new ChunkBookmarkData;
    VoxelData * __restrict voxel_data0 = new VoxelData;
    VoxelData * __restrict voxel_data1 = new VoxelData;

    ChunkBookmarkData * __restrict all_chunk_data = chunk_data0;
    VoxelData * __restrict all_voxel_data = voxel_data0;

    VoxelRenderData * __restrict voxel_render_data = new VoxelRenderData;

    f32 frame_timer = 0.0f;
    f32 last_time = 0.0f;
    const static f32 TIME_STEP = 0.016f;
    b32 running = true;
    while(running)
    {
        f64 current_time = glfwGetTime();
        frame_timer += (f32)current_time - last_time;
        last_time = current_time;

        if(frame_timer >= TIME_STEP)
        {
            frame_timer -= TIME_STEP;
            glfwPollEvents();

            // Begin frame
            {
                glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                ImGui_ImplOpenGL3_NewFrame();
                ImGui_ImplGlfw_NewFrame();
                ImGui::NewFrame();
            }

            running = !glfwGetKey(graphics_state->window, GLFW_KEY_ESCAPE) && !glfwWindowShouldClose(graphics_state->window);

            {
                ImGui::Text("camera pos (%f, %f, %f)", graphics_state->cam.pos.x, graphics_state->cam.pos.y, graphics_state->cam.pos.z);
                ImGui::Text("camera x rot %f", graphics_state->cam.rot_x);
                ImGui::Text("camera y rot %f", graphics_state->cam.rot_y);
                ImGui::DragFloat("camera near_plane_dist", &graphics_state->cam.near_plane_dist);
                ImGui::DragFloat("camera view_dist", &graphics_state->cam.view_dist);
                ImGui::Checkbox("debug cam", &graphics_state->use_debug_cam);

                Camera *current_cam;
                if(graphics_state->use_debug_cam) current_cam = &graphics_state->debug_cam;
                else current_cam = &graphics_state->cam;

                static f32 last_mouse_x, last_mouse_y;
                f64 xpos, ypos;
                glfwGetCursorPos(graphics_state->window, &xpos, &ypos);
                if(glfwGetKey(graphics_state->window, GLFW_KEY_SPACE))
                {
                    v2 mouse_delta = v2((f32)xpos - last_mouse_x, (f32)ypos - last_mouse_y);
                    current_cam->rot_y -= mouse_delta.x * TIME_STEP * 0.2f;
                    current_cam->rot_x -= mouse_delta.y * TIME_STEP * 0.2f;
                }
                last_mouse_x = (f32)xpos;
                last_mouse_y = (f32)ypos;

                mat4 cam_xform = view_m_world(current_cam);
                v3 cam_i = v3(cam_xform[0][0], cam_xform[0][1], cam_xform[0][2]);
                v3 cam_k = v3(cam_xform[2][0], cam_xform[2][1], cam_xform[2][2]);

                v3 camera_vel = v3();
                camera_vel += -cam_i * glfwGetKey(graphics_state->window, GLFW_KEY_A);
                camera_vel +=  cam_i * glfwGetKey(graphics_state->window, GLFW_KEY_D);
                camera_vel +=  cam_k * glfwGetKey(graphics_state->window, GLFW_KEY_S);
                camera_vel += -cam_k * glfwGetKey(graphics_state->window, GLFW_KEY_W);
                static f32 camera_speed = 100.0f;
                ImGui::DragFloat("camera speed", &camera_speed);
                current_cam->pos += normalize(camera_vel) * TIME_STEP * camera_speed;
            }


            s32 camera_x = (s32)graphics_state->cam.pos.x;
            s32 camera_y = (s32)graphics_state->cam.pos.y;
            s32 camera_z = (s32)graphics_state->cam.pos.z;
            s32 camera_chunk_x = camera_x >> CHUNK_POW;
            s32 camera_chunk_y = camera_y >> CHUNK_POW;
            s32 camera_chunk_z = camera_z >> CHUNK_POW;

            s32 region_chunks_bl_x = camera_chunk_x - LOADED_REGION_CHUNKS_DIM/2;
            s32 region_chunks_bl_y = camera_chunk_y - LOADED_REGION_CHUNKS_DIM/2;
            s32 region_chunks_bl_z = camera_chunk_z - LOADED_REGION_CHUNKS_DIM/2;
            s32 region_chunks_tr_x = camera_chunk_x + LOADED_REGION_CHUNKS_DIM/2;
            s32 region_chunks_tr_y = camera_chunk_y + LOADED_REGION_CHUNKS_DIM/2;
            s32 region_chunks_tr_z = camera_chunk_z + LOADED_REGION_CHUNKS_DIM/2;

            //static s32 last_camera_chunk_x = 0xFFFFFFFF;
            //static s32 last_camera_chunk_y = 0xFFFFFFFF;
            //static s32 last_camera_chunk_z = 0xFFFFFFFF;
            static s32 last_camera_chunk_x = camera_chunk_x;
            static s32 last_camera_chunk_y = camera_chunk_y;
            static s32 last_camera_chunk_z = camera_chunk_z;

            static s64 last_gen_time = 0;

            // Initial terrain generation around player
            static bool do_gen = true;
            if(do_gen)
            {
                all_chunk_data->num = 0;
                for(s32 x = region_chunks_bl_x; x < region_chunks_tr_x; x++)
                {
                    for(s32 y = region_chunks_bl_y; y < region_chunks_tr_y; y++)
                    {
                        for(s32 z = region_chunks_bl_z; z < region_chunks_tr_z; z++)
                        {
                            gen_chunk(
                                all_chunk_data,
                                all_voxel_data,
                                x, y, z
                            );
                        }
                    }
                }
                // HACK
                all_voxel_data->num = all_chunk_data->bookmark[all_chunk_data->num];
                do_gen = false;
            }

            ChunkBookmarkData* src_chunk_data = nullptr;
            VoxelData* src_voxel_data = nullptr;
            ChunkBookmarkData* dst_chunk_data = nullptr;
            VoxelData* dst_voxel_data = nullptr;
            if(all_chunk_data == chunk_data0)
            {
                src_chunk_data = chunk_data0;
                src_voxel_data = voxel_data0;
                dst_chunk_data = chunk_data1;
                dst_voxel_data = voxel_data1;
            }
            else
            {
                src_chunk_data = chunk_data1;
                src_voxel_data = voxel_data1;
                dst_chunk_data = chunk_data0;
                dst_voxel_data = voxel_data0;
            }
            bool did_generation = maybe_gen_new_terrain(
                dst_chunk_data,
                dst_voxel_data,
                src_chunk_data,
                src_voxel_data,

                last_camera_chunk_x,
                last_camera_chunk_y,
                last_camera_chunk_z,
                camera_chunk_x,
                camera_chunk_y,
                camera_chunk_z
            );
            if(did_generation)
            {
                if(all_chunk_data == chunk_data0)
                {
                    all_chunk_data = chunk_data1;
                    all_voxel_data = voxel_data1;
                }
                else
                {
                    all_chunk_data = chunk_data0;
                    all_voxel_data = voxel_data0;
                }
            }

            last_camera_chunk_x = camera_chunk_x;
            last_camera_chunk_y = camera_chunk_y;
            last_camera_chunk_z = camera_chunk_z;
            ImGui::Text("Last gen time: %ims", last_gen_time);


            // Prep render data
            {
                TIME_SCOPE("Prep render data");

                ImGui::Text("prep num_voxels %i", all_voxel_data->num);

                camera_cull(
                    &voxel_render_data->num,
                    voxel_render_data->pos,
                    voxel_render_data->color,
                    all_voxel_data->num,
                    all_voxel_data->pos,
                    all_voxel_data->color,
                    &graphics_state->cam
                );
            }

            // Draw
            {
                TIME_SCOPE("Draw");

                u32 num_batches_drawn = 0;
                for(u32 batch_i = 0; batch_i < voxel_render_data->num; batch_i += VoxelRenderData::BATCH_SIZE)
                {
                    f32 *batch_x = voxel_render_data->pos + MAX_VOXELS*0 + batch_i;
                    f32 *batch_y = voxel_render_data->pos + MAX_VOXELS*1 + batch_i;
                    f32 *batch_z = voxel_render_data->pos + MAX_VOXELS*2 + batch_i;
                    u32 *batch_color = voxel_render_data->color + batch_i;
                    u32 batch_size = min(voxel_render_data->num - batch_i, VoxelRenderData::BATCH_SIZE);

                    assert(VoxelRenderData::BATCH_SIZE <= SHADER_BUFFER_WIDTH);
                    glUseProgram(graphics_state->batch_voxel_shader_program);
                    check_gl_errors("use program");

                    Camera *current_cam;
                    if(graphics_state->use_debug_cam) current_cam = &graphics_state->debug_cam;
                    else current_cam = &graphics_state->cam;

                    {
                        mat4 m = view_m_world(current_cam);
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "m_view");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }
                    {
                        mat4 m = clip_m_view(current_cam);
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "m_proj");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }

                    glBindVertexArray(graphics_state->batch_voxel_vao);
                    glBindBuffer(GL_ARRAY_BUFFER, graphics_state->batch_voxel_vbo);
                    u32 offset_x = 0 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_y = 1 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_z = 2 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_color = 3 * SHADER_BUFFER_WIDTH*sizeof(u32);
                    u32 size = sizeof(f32) * batch_size;
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, graphics_state->batch_voxel_ssbo);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_x, size, batch_x);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_y, size, batch_y);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_z, size, batch_z);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_color, size, batch_color);
                    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, graphics_state->batch_voxel_ssbo);

                    glDrawElementsInstanced(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0, batch_size);
                    check_gl_errors("draw");

                    num_batches_drawn++;
                }
                ImGui::Text("num_voxels %i", voxel_render_data->num);
                ImGui::Text("batches %i", num_batches_drawn);

            }

            {
                draw_line(v3(), v3(1.0f, 0.0f, 0.0f), v3(1.0f, 0.0f, 0.0f));
                draw_line(v3(), v3(0.0f, 1.0f, 0.0f), v3(0.0f, 1.0f, 0.0f));
                draw_line(v3(), v3(0.0f, 0.0f, 1.0f), v3(0.0f, 0.0f, 1.0f));

                mat4 cam_xform = view_m_world(&graphics_state->cam);
                v3 cam_p = graphics_state->cam.pos;
                v3 cam_i = v3(cam_xform[0][0], cam_xform[0][1], cam_xform[0][2]);
                v3 cam_j = v3(cam_xform[1][0], cam_xform[1][1], cam_xform[1][2]);
                v3 cam_k = v3(cam_xform[2][0], cam_xform[2][1], cam_xform[2][2]);
                draw_line(cam_p, cam_p + cam_i, v3(1.0f, 0.0f, 0.0f));
                draw_line(cam_p, cam_p + cam_j, v3(0.0f, 1.0f, 0.0f));
                draw_line(cam_p, cam_p + cam_k, v3(0.0f, 0.0f, 1.0f));

                f32 tan_fov = tan(graphics_state->cam.vfov*0.5f);
                f32 near_hh = tan_fov * graphics_state->cam.near_plane_dist;
                f32 near_hw = near_hh * graphics_state->cam.ar;
                f32 far_hh = tan_fov * graphics_state->cam.view_dist;
                f32 far_hw = far_hh * graphics_state->cam.ar;

                v3 p_near = cam_p + -cam_k * graphics_state->cam.near_plane_dist;
                v3 p_far = cam_p + -cam_k * graphics_state->cam.view_dist;

                v3 frustum_v[] =
                {
                    p_near - cam_i*near_hw - cam_j*near_hh, // bl 
                    p_near + cam_i*near_hw - cam_j*near_hh, // br 
                    p_near + cam_i*near_hw + cam_j*near_hh, // tr 
                    p_near - cam_i*near_hw + cam_j*near_hh, // tl 

                    p_far - cam_i*far_hw - cam_j*far_hh, // bl 
                    p_far + cam_i*far_hw - cam_j*far_hh, // br 
                    p_far + cam_i*far_hw + cam_j*far_hh, // tr 
                    p_far - cam_i*far_hw + cam_j*far_hh, // tl 
                };

                draw_line(frustum_v[0], frustum_v[1], v3(1.0f, 0.0f, 0.0f));
                draw_line(frustum_v[1], frustum_v[2], v3(1.0f, 0.0f, 0.0f));
                draw_line(frustum_v[2], frustum_v[3], v3(1.0f, 0.0f, 0.0f));
                draw_line(frustum_v[3], frustum_v[0], v3(1.0f, 0.0f, 0.0f));

                draw_line(frustum_v[4], frustum_v[5], v3(0.0f, 1.0f, 0.0f));
                draw_line(frustum_v[5], frustum_v[6], v3(0.0f, 1.0f, 0.0f));
                draw_line(frustum_v[6], frustum_v[7], v3(0.0f, 1.0f, 0.0f));
                draw_line(frustum_v[7], frustum_v[4], v3(0.0f, 1.0f, 0.0f));

                // Right
                draw_line(frustum_v[1], frustum_v[5], v3(0.0f, 0.0f, 1.0f));
                draw_line(frustum_v[5], frustum_v[6], v3(0.0f, 0.0f, 1.0f));
                draw_line(frustum_v[6], frustum_v[2], v3(0.0f, 0.0f, 1.0f));
                draw_line(frustum_v[2], frustum_v[1], v3(0.0f, 0.0f, 1.0f));

                // Top
                draw_line(frustum_v[2], frustum_v[6], v3(1.0f, 0.0f, 1.0f));
                draw_line(frustum_v[6], frustum_v[7], v3(1.0f, 0.0f, 1.0f));
                draw_line(frustum_v[7], frustum_v[3], v3(1.0f, 0.0f, 1.0f));
                draw_line(frustum_v[3], frustum_v[2], v3(1.0f, 0.0f, 1.0f));

                // Left
                draw_line(frustum_v[3], frustum_v[7], v3(0.0f, 1.0f, 1.0f));
                draw_line(frustum_v[7], frustum_v[4], v3(0.0f, 1.0f, 1.0f));
                draw_line(frustum_v[4], frustum_v[0], v3(0.0f, 1.0f, 1.0f));
                draw_line(frustum_v[0], frustum_v[3], v3(0.0f, 1.0f, 1.0f));

                // Bottom
                draw_line(frustum_v[0], frustum_v[4], v3(1.0f, 1.0f, 0.0f));
                draw_line(frustum_v[4], frustum_v[5], v3(1.0f, 1.0f, 0.0f));
                draw_line(frustum_v[5], frustum_v[1], v3(1.0f, 1.0f, 0.0f));
                draw_line(frustum_v[1], frustum_v[0], v3(1.0f, 1.0f, 0.0f));
            }

            static bool show_chunk_lines = true;
            ImGui::Checkbox("show_chunk_lines", &show_chunk_lines);
            if(show_chunk_lines)
            {
                for(s32 x = region_chunks_bl_x; x < region_chunks_tr_x; ++x)
                {
                    for(s32 y = region_chunks_bl_y; y < region_chunks_tr_y; ++y)
                    {
                        for(s32 z = region_chunks_bl_z; z < region_chunks_tr_z; ++z)
                        {
                            draw_line(v3(f32(x*CHUNK_DIM),       f32(y*CHUNK_DIM), f32(z*CHUNK_DIM)),
                                      v3(f32((x + 1)*CHUNK_DIM), f32(y*CHUNK_DIM), f32(z*CHUNK_DIM)),
                                      v3(1.0f, 0.0f, 0.0f));
                            draw_line(v3(f32(x*CHUNK_DIM), f32(y*CHUNK_DIM),       f32(z*CHUNK_DIM)),
                                      v3(f32(x*CHUNK_DIM), f32((y + 1)*CHUNK_DIM), f32(z*CHUNK_DIM)),
                                      v3(0.0f, 1.0f, 0.0f));
                            draw_line(v3(f32(x*CHUNK_DIM), f32(y*CHUNK_DIM), f32(z*CHUNK_DIM)),
                                      v3(f32(x*CHUNK_DIM), f32(y*CHUNK_DIM), f32((z + 1)*CHUNK_DIM)),
                                      v3(0.0f, 0.0f, 1.0f));
                        }
                    }
                }
            }

            {
                f64 frame_ms = (glfwGetTime() - current_time) * 1000.0;
                static f64 frame_times[32] = {};
                static u32 num_frame_times = 0;
                frame_times[num_frame_times++] = frame_ms;
                if(num_frame_times >= 32) num_frame_times = 0;
                f64 frame_time_sum = 0.0;
                for(u32 i = 0; i < 32; i++) frame_time_sum += frame_times[i];
                ImGui::Text("frame time: %f ms", frame_time_sum / 32.0);

                ImGui::Text(graphics_state->debug_info_log);
            }

            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(graphics_state->window);
        }
    }

    {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();

        glfwTerminate();
    }

    return 0;
}

