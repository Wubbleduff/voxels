
#include <windows.h>
#include <stdio.h> // fopen
#include <immintrin.h>
#include <atomic>

// std sort
#include <algorithm>

#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "common.h"
#include "game_math.h"

// TODO(mfritz) delete
#include <algorithm>

#pragma warning(disable:4996)

// Globals
// NOTE: This value must match the shader.
static constexpr u32 SHADER_BUFFER_WIDTH = 10*1024;
static constexpr u64 MAX_VOXELS = 50*1024*1024;
static constexpr u32 NUM_TOTAL_THREADS = 16;
static constexpr u32 MAX_CHUNKS = 1 << 17;
static constexpr u32 MAX_LOD = 11;
//static constexpr u32 MAX_CHUNKS_GEN_PER_FRAME = MAX_CHUNKS;
static constexpr u32 MAX_CHUNKS_GEN_PER_FRAME = 32;
static struct GraphicsState *G_graphics_state = nullptr;
static struct InputState *G_input_state = nullptr;
static FILE* G_log_file = nullptr;

#define ASSERT(exp, msg)                        \
{                                           \
    if(!(exp))                              \
    {                                       \
        fprintf(G_log_file, "%s:%i @ %s @ %s", __FILE__, __LINE__, #exp, msg); \
        fflush(G_log_file);                 \
        DebugBreak();                       \
    }                                       \
}

struct StringBuf
{
    static constexpr u32 SIZE = 32;
    char s[SIZE];
};
static StringBuf make_string_buf(const char* s)
{
    ASSERT(strlen(s) < 32, "String too long for StringBuf");
    StringBuf r;
    strncpy(r.s, s, StringBuf::SIZE);
    return r;
}

////////////////////////////////////////////////////////////////////////////////
// Graphics
////////////////////////////////////////////////////////////////////////////////
struct VoxelRenderData
{
    u32 num = 0;
    // | xxxx yyyy zzz |
    f32 pos[MAX_VOXELS*3];
    f32 scale[MAX_VOXELS];
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
static constexpr f32 SQRT3 = 1.73205080f;
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
    
    u32 imgui_debug_texture_width;
    u32 imgui_debug_texture_height;
    GLuint imgui_debug_texture;
    u32* imgui_debug_texture_data;
    
    Camera cam;
    Camera debug_cam;
    bool use_debug_cam;
    
    //static constexpr f32 voxel_vertices[32] =
    static constexpr f32 voxel_vertices[] =
    {
        0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, // X
        0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, // Y
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, // Z
        -SQRT3,  SQRT3, -SQRT3,  SQRT3, -SQRT3,  SQRT3, -SQRT3,  SQRT3, 
        -SQRT3, -SQRT3,  SQRT3,  SQRT3, -SQRT3, -SQRT3,  SQRT3,  SQRT3, 
        -SQRT3, -SQRT3, -SQRT3, -SQRT3,  SQRT3,  SQRT3,  SQRT3,  SQRT3, 
    };

    static constexpr u32 voxel_indices[] = 
    {
        0, 4, 2, // -X
        2, 4, 6, // -X

        3, 5, 1, // +X
        3, 7, 5, // +X

        0, 1, 4, // -Y
        4, 1, 5, // -Y

        2, 6, 7, // +Y
        3, 2, 7, // +Y

        1, 0, 2, // -Z
        2, 3, 1, // -Z

        4, 5, 6, // +Z
        6, 5, 7, // +Z
    };
};
struct InputState
{
    s32 mouse_screen_x;
    s32 mouse_screen_y;
};

// TODO speed
mat4 view_m_world(const Camera* cam)
{
    mat4 y_rot = make_y_axis_rotation_matrix(-cam->rot_y);
    mat4 x_rot = make_x_axis_rotation_matrix(-cam->rot_x);
    mat4 result = x_rot * y_rot * make_translation_matrix(-cam->pos);
    return result;
}

mat4 clip_m_view(const Camera* cam)
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

void glfw_error_fn(s32 error, const char *message)
{
    (void)error;
    (void)message;
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

////////////////////////////////////////////////////////////////////////////////
// Voxel data
////////////////////////////////////////////////////////////////////////////////

// Indicates voxel is empty for dense arrays.
// Voxel IDs are valid for values [0, EMPTY_VOXEL)
static constexpr u32 EMPTY_VOXEL = u32(-1);

struct Chunk
{
    v3i bl_pos;
    u32 lod;
    u32 num;
    // VLA | xxxx yyyy zzzz cccc |
    static constexpr u32 INTS_PER_VOXEL = 4;
    s32 voxels[1];
};

struct PackedVoxels
{
    u32 num;
    s32 pos[MAX_VOXELS*3];
    // TODO(mfritz) Does this need to be per voxel? Can we just group by chunks?
    s32 scale[MAX_VOXELS];
    u32 voxel_id[MAX_VOXELS];
};

// Sparse voxel array.
struct TreeBuffer
{
    static constexpr u32 CAP = 32768;
    u32 num;
    s32 pos[CAP*3]; // xxx yyy zzz
    u32 voxel_id[CAP];
};
struct Tree
{
    static constexpr u32 MAX_NODES = 8;
    static constexpr u32 MAX_BRANCHES = 8;
    
    u32 num_nodes;
    v3 nodes_pos[MAX_NODES];
    f32 nodes_dist[MAX_NODES];
    
    u32 num_branches;
    Tree* branches[8];
    
    v3 get_point_at_height(f32 h)
    {
        ASSERT(num_nodes > 0, "Can't find point on empty tree.");
        if(num_nodes == 1)
        {
            ASSERT(h == 0.0f, "Tree with only 1 node must have h of 0.");
            return nodes_pos[0];
        }
        
        u32 node_idx = 0;
        while(h > nodes_dist[node_idx + 1])
        {
            node_idx++;
            ASSERT(node_idx < num_nodes, "Finding point on tree given h is too tall.");
        }
        
        f32 t = (h - nodes_dist[node_idx]) / (nodes_dist[node_idx + 1] - nodes_dist[node_idx]);
        return lerp(nodes_pos[node_idx], nodes_pos[node_idx + 1], t);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Profiling
////////////////////////////////////////////////////////////////////////////////

struct ProfileData
{
    static constexpr u32 MAX_TIMES = 64;
    u32 num;
    s64 time[MAX_TIMES];
    StringBuf name[MAX_TIMES];
};
static ProfileData G_profile_data;

void record_profile_data(const char* name, s64 t)
{
    u32& num = G_profile_data.num;
    G_profile_data.time[num] = t;
    G_profile_data.name[num] = make_string_buf(name);
    num++;
}

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
        ms.QuadPart = (end.QuadPart - start.QuadPart) * 1'000'000'000;
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
        
        LARGE_INTEGER ns;
        ns.QuadPart = (end.QuadPart - start.QuadPart) * 1'000'000'000;
        ns.QuadPart /= freq.QuadPart;
        
        //ImGui::Text("TIME - %s: %ius", m_name, ms.QuadPart);
        record_profile_data(m_name, ns.QuadPart / 1'000);
    }
    const char *m_name;
    LARGE_INTEGER start;
};
#define TIME_SCOPE(name) ScopeTimer _time_scope = ScopeTimer(name)

////////////////////////////////////////////////////////////////////////////////
// EXPERIMENTAL
////////////////////////////////////////////////////////////////////////////////
struct CubeIterator
{
    s32 idx = 0;
    s32 x = 0;
    s32 y = 0;
    s32 z = 0;
};
#define ITER_CUBE(N)                                             \
for(CubeIterator it;                                         \
it.idx < (N)*(N)*(N);                                    \
it.idx++,                                                \
it.z = it.idx / ((N)*(N)),                               \
it.y = (it.idx / (N)) % (N),                             \
it.x = it.idx % (N) )                                    


#define KB(n) ((n) / 1024)
#define MB(n) ((n) / 1024 / 1024)


////////////////////////////////////////////////////////////////////////////////
// Job System
////////////////////////////////////////////////////////////////////////////////
struct ThreadState
{
    Chunk* gen_chunk_scratch = nullptr;
};

static ThreadState* main_thread_state;

enum class JobId : u64
{
    load_terrain,
    
    num_jobs
};
typedef void (*WorkerFn)(ThreadState*, void**, const void*);
constexpr u32 MAX_WORK_ITEMS = 1024*1024;
struct WorkQueue
{
    std::atomic<u32> in_begin = 0;
    std::atomic<u32> in_end = 0;
    WorkerFn in_fn[MAX_WORK_ITEMS];
    const void* in_data[MAX_WORK_ITEMS];
    
    std::atomic<u32> out_begin = 0;
    std::atomic<u32> out_end = 0;
    void** out_data[MAX_WORK_ITEMS];
    
    std::atomic<u32> pending_work_count = 0;
};
static WorkQueue G_work_queue_load_terrain;
WorkQueue* get_queue_from_job_id(JobId job_id)
{
    WorkQueue* q = nullptr;
    switch(job_id)
    {
        case JobId::load_terrain:
        q = &G_work_queue_load_terrain;
        break;
        default:
        assert(false);
        break;
    }
    return q;
}
DWORD WINAPI worker_main(LPVOID lpParam)
{
    ThreadState* thread_state = reinterpret_cast<ThreadState*>(lpParam);
    
    for(u64 job_id = 0; job_id < u64(JobId::num_jobs); job_id++)
    {
        // TODO(mfritz) This can just be stored in a flat array rather than looking up each time.
        WorkQueue* q = get_queue_from_job_id(JobId(job_id));
        while(true)
        {
            u32 q_begin = q->in_begin.load();
            const u32 q_end = q->in_end.load();
            if(q_begin == q_end)
            {
                _mm_pause();
                continue;
            }
            if(q->in_begin.compare_exchange_weak(q_begin, q_begin + 1))
            {
                q->in_fn[q_begin](thread_state, q->out_data[q_begin], q->in_data[q_begin]);
                q->pending_work_count--;
            }
        }
    }
    
    return 0;
}
// Must only be called from the main thread.
void add_work(void** out_data, const JobId job_id, const WorkerFn fn, const void* data)
{
    WorkQueue* q = get_queue_from_job_id(job_id);
    u32 q_end = q->in_end.load();
    q->in_fn[q_end] = fn;
    q->in_data[q_end] = data;
    q->out_data[q_end] = out_data;
    q->in_end++;
    q->pending_work_count++;
}
// Must only be called from the main thread.
void wait_for_job(const JobId job_id)
{
    WorkQueue* q = get_queue_from_job_id(job_id);
    while(true)
    {
        _mm_pause();
        u32 pending_work_count = q->pending_work_count.load();
        if(pending_work_count == 0)
        {
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Debug draw
////////////////////////////////////////////////////////////////////////////////

static void draw_line(v3 a, v3 b, v3 c)
{
    glUseProgram(G_graphics_state->debug_line_shader_program);
    check_gl_errors("use program");
    
    Camera *current_cam;
    if(G_graphics_state->use_debug_cam) current_cam = &G_graphics_state->debug_cam;
    else current_cam = &G_graphics_state->cam;
    
    mat4 m = clip_m_view(current_cam) * view_m_world(current_cam);
    GLint loc = glGetUniformLocation(G_graphics_state->debug_line_shader_program, "mvp");
    glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
    if(loc == -1) assert(false);
    
    loc = glGetUniformLocation(G_graphics_state->debug_line_shader_program, "color");
    glUniform3f(loc, c.r, c.g, c.b);
    if(loc == -1) assert(false);
    
    glBindVertexArray(G_graphics_state->debug_line_vao);
    glBindBuffer(GL_ARRAY_BUFFER, G_graphics_state->debug_line_vbo);
    
    f32 v[] =
    {
        a.x, a.y, a.z,
        b.x, b.y, b.z
    };
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(v), v);
    
    glDrawArrays(GL_LINES, 0, 2);
    check_gl_errors("draw");
}

static void debug_draw_frustum()
{
    TIME_SCOPE("debug draw frustum");
    
    draw_line(v3(), v3(1.0f, 0.0f, 0.0f), v3(1.0f, 0.0f, 0.0f));
    draw_line(v3(), v3(0.0f, 1.0f, 0.0f), v3(0.0f, 1.0f, 0.0f));
    draw_line(v3(), v3(0.0f, 0.0f, 1.0f), v3(0.0f, 0.0f, 1.0f));
    
    mat4 cam_xform = view_m_world(&G_graphics_state->cam);
    v3 cam_p = G_graphics_state->cam.pos;
    v3 cam_i = v3(cam_xform[0][0], cam_xform[0][1], cam_xform[0][2]);
    v3 cam_j = v3(cam_xform[1][0], cam_xform[1][1], cam_xform[1][2]);
    v3 cam_k = v3(cam_xform[2][0], cam_xform[2][1], cam_xform[2][2]);
    draw_line(cam_p, cam_p + cam_i, v3(1.0f, 0.0f, 0.0f));
    draw_line(cam_p, cam_p + cam_j, v3(0.0f, 1.0f, 0.0f));
    draw_line(cam_p, cam_p + cam_k, v3(0.0f, 0.0f, 1.0f));
    
    f32 tan_fov = tan(G_graphics_state->cam.vfov*0.5f);
    f32 near_hh = tan_fov * G_graphics_state->cam.near_plane_dist;
    f32 near_hw = near_hh * G_graphics_state->cam.ar;
    f32 far_hh = tan_fov * G_graphics_state->cam.view_dist;
    f32 far_hw = far_hh * G_graphics_state->cam.ar;
    
    v3 p_near = cam_p + -cam_k * G_graphics_state->cam.near_plane_dist;
    v3 p_far = cam_p + -cam_k * G_graphics_state->cam.view_dist;
    
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

////////////////////////////////////////////////////////////////////////////////
// Code
////////////////////////////////////////////////////////////////////////////////

// Outputs plan normals for right, top, left, bottom, near, far in that order.
static void get_frustum_planes(v3* out_normals, v3* out_points)
{
    mat4 cam_xform = view_m_world(&G_graphics_state->cam);
    v3 cam_p = G_graphics_state->cam.pos;
    v3 cam_i = v3(cam_xform[0][0], cam_xform[0][1], cam_xform[0][2]);
    v3 cam_j = v3(cam_xform[1][0], cam_xform[1][1], cam_xform[1][2]);
    v3 cam_k = v3(cam_xform[2][0], cam_xform[2][1], cam_xform[2][2]);
    
    f32 tan_fov = tan(deg_to_rad(G_graphics_state->cam.vfov)*0.5f);
    f32 near_hh = tan_fov * G_graphics_state->cam.near_plane_dist;
    f32 near_hw = near_hh * G_graphics_state->cam.ar;
    f32 far_hh = tan_fov * G_graphics_state->cam.view_dist;
    f32 far_hw = far_hh * G_graphics_state->cam.ar;
    
    v3 p_near = cam_p + -cam_k * G_graphics_state->cam.near_plane_dist;
    v3 p_far = cam_p + -cam_k * G_graphics_state->cam.view_dist;
    
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
    // TODO(mfritz) Think harder about this...
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
                 u32* out_num_voxels,
                 f32* out_voxel_pos,
                 f32* out_voxel_scale,
                 u32* out_voxel_id,
                 
                 const u32 in_num_voxels,
                 const s32* in_voxel_pos,
                 const s32* in_voxel_scale,
                 const u32* in_voxel_id)
{
    // Check distance from voxel centroid to each plane and compare against 1.0f.
    // This is effectively the same as a sphere-plane distance check.
    // If P is a point on the plane, Q is the voxel centroid, and n is the plane normal, dist check is:
    // (Q - P) * n < sqrt(0.5^2 + 0.5^2) ->
    // Q*n - P*n - sqrt(0.5^2 + 0.5^2)
    // We can precompute p_dot_n = P*n + sqrt(0.5^2 + 0.5^2) outside the per-voxel loop.
    // Inside the loop, we just need to compute Q*n - p_dot_n. The resulting dist will
    // have the sign bit set if < 0. Then we can use movemask without needing an extra
    // cmp.
    
    v3 plane_normals[6];
    v3 plane_points[6];
    get_frustum_planes(plane_normals, plane_points);
    
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
        __m256 vs = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(in_voxel_scale + i)));
        __m256i vc = _mm256_loadu_si256((__m256i*)(in_voxel_id + i));
        
        __m256 mask = _mm256_cvtepi32_ps(_mm256_set1_epi8(-1));
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
        _mm256_storeu_ps(out_voxel_scale + num_voxels, left_pack(vs, pack_index));
        _mm256_storeu_si256((__m256i*)(out_voxel_id + num_voxels), left_pack(vc, pack_index));
        num_voxels += _mm_popcnt_u32(pack_index);
    }
    *out_num_voxels = num_voxels;
}

static void reshape(GLFWwindow *window, s32 width, s32 height)
{
    (void)window;
    glViewport(0, 0, (GLint)width, (GLint)height);
    G_graphics_state->cam.ar = (f32)width / height;
    G_graphics_state->debug_cam.ar = (f32)width / height;
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
    
    u64 current_tag_length = 0;
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
    const u32 info_log_size = sizeof(G_graphics_state->debug_info_log);
    
    u32 vert_shader;
    vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vert_source, NULL);
    glCompileShader(vert_shader);
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(vert_shader, info_log_size, NULL, G_graphics_state->debug_info_log);
        return GLuint(-1);
    }
    
    u32 frag_shader;
    frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(frag_shader, 1, &frag_source, NULL);
    glCompileShader(frag_shader);
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(frag_shader, info_log_size, NULL, G_graphics_state->debug_info_log);
        return GLuint(-1);
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
        glGetProgramInfoLog(program, info_log_size, NULL, G_graphics_state->debug_info_log);
        return GLuint(-1);
    }
    
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader); 
    check_gl_errors("deleting shaders");
    
    return program;
}

GLuint make_shader_from_file(const char *shader_path)
{
    char *vert_source = nullptr;
    char *geom_source = nullptr;
    char *frag_source = nullptr;
    bool read = read_shader_file(shader_path, &vert_source, &geom_source, &frag_source);
    if(!read)
    {
        fprintf(G_log_file, "Could not open shader file %s\n", shader_path);
        return GLuint(-1);
    }
    
    return make_shader_from_string(vert_source, frag_source);
}

bool point_cube_intersect(v3i p, v3i bl, s32 dim)
{
    bool result =
        !(p.x() < bl.x() || p.x() >= bl.x() + dim ||
          p.y() < bl.y() || p.y() >= bl.y() + dim ||
          p.z() < bl.z() || p.z() >= bl.z() + dim);
    return result;
}

bool cube_cube_intersect(v3i bl0, v3i bl1, s32 dim)
{
    bool result =
        !(bl0.x() + dim < bl1.x() || bl0.x() >= bl1.x() + dim ||
          bl0.y() + dim < bl1.y() || bl0.y() >= bl1.y() + dim ||
          bl0.z() + dim < bl1.z() || bl0.z() >= bl1.z() + dim);
    return result;
}

static v3i point_to_voxel(v3 p)
{
    __m128 q = _mm_setr_ps(p.x, p.y, p.z, 0.0f);
    q = _mm_round_ps(q, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
    v3i r(_mm_cvtps_epi32(q));
    return r;
}

__m256 sample_terrain(__m256 x, __m256 y, __m256 z)
{
    __m256 total_n = _mm256_set1_ps(0.0f);

    {
        __m256 sx = _mm256_mul_ps(x, _mm256_set1_ps(0.006f));
        __m256 sy = _mm256_mul_ps(y, _mm256_set1_ps(0.040f));
        __m256 sz = _mm256_mul_ps(z, _mm256_set1_ps(0.006f));
        __m256 n = pnoise8(sx, sy, sz);

        n = _mm256_mul_ps(n, _mm256_set1_ps(0.7f));

        constexpr f32 squashing_factor = 0.25f;
        n = _mm256_sub_ps(n, _mm256_mul_ps(sy, _mm256_set1_ps(squashing_factor)));

        total_n = _mm256_add_ps(total_n, n);
    }

    {
        __m256 sx = _mm256_mul_ps(x, _mm256_set1_ps(0.0005f));
        __m256 sy = _mm256_mul_ps(y, _mm256_set1_ps(0.0004f));
        __m256 sz = _mm256_mul_ps(z, _mm256_set1_ps(0.0005f));
        __m256 n = pnoise8(sx, sy, sz);

        n = _mm256_mul_ps(n, n);
        n = _mm256_mul_ps(n, _mm256_set1_ps(200.0f));

        total_n = _mm256_add_ps(total_n, n);
    }

    return total_n;
}


Chunk* allocate_chunk(u32 cap_voxels)
{
    u64 bytes = sizeof(Chunk) + cap_voxels * sizeof(s32) * Chunk::INTS_PER_VOXEL;
    Chunk* result = reinterpret_cast<Chunk*>(_mm_malloc(bytes, 4096));
    return result;
}

void deallocate_chunk(Chunk** chunk)
{
    _mm_free(*chunk);
    *chunk = nullptr;
}


static constexpr u32 MIN_CHUNK_DIM = 32;
// Each chunk will be 32x32x32 voxels. The scale of each voxel in the chunk
// is determined by its LOD.
static constexpr u32 VOXELS_CHUNK_DIM = 32;
static constexpr u32 CHUNK_MAX_VOXELS = VOXELS_CHUNK_DIM*VOXELS_CHUNK_DIM*VOXELS_CHUNK_DIM;
static u32 get_chunk_lod_from_distance(u32 distance)
{
    //return distance / 1024;
    //return distance / 256;
    //return u32(sqrt(double(distance)));

    if(distance < 1 << 7) return 0;
    else if(distance < 1 << 8)  return 0;
    else if(distance < 1 << 9)  return 1;
    else if(distance < 1 << 10) return 2;
    else if(distance < 1 << 11) return 3;
    else if(distance < 1 << 12) return 4;
    else if(distance < 1 << 13) return 5;
    else if(distance < 1 << 14) return 6;
    else return 7;
}
static u32 get_chunk_dim_from_lod(u32 lod)
{
    return (1 << lod) * MIN_CHUNK_DIM;
}
static u32 get_lod_scale(u32 lod)
{
    return 1 << lod;
}
static u32 get_chunk_target_lod(v3i chunk_bl_pos, u32 cur_chunk_lod, v3i player_pos)
{
    // Given a chunk and the player's position, calculate what thet target LOD should be for that chunk.
    // A chunk's LOD is only a function of the chunk's position, size, and its distance from the player.
    // The distance calc is defined as the distance from the bounding sphere's surface to the player point
    // (clamped to 0 for negative dist).
    const u32 chunk_pos_mask = ~((1U << 5U) - 1U);
    static_assert(MIN_CHUNK_DIM == 1U << 5U);
    const v3i player_chunk_pos = _mm_and_si128(player_pos.v, v3i(chunk_pos_mask, chunk_pos_mask, chunk_pos_mask).v);
    u32 chunk_hdim = get_chunk_dim_from_lod(cur_chunk_lod) / 2;
    v3i chunk_center = chunk_bl_pos + v3i(chunk_hdim, chunk_hdim, chunk_hdim);
    u32 dist_to_player = u32(dot(chunk_center - player_chunk_pos, chunk_center - player_chunk_pos));
    dist_to_player = u32(sqrtf(float(dist_to_player)));
    u32 chunk_radius = u32(double(chunk_hdim*2) * SQRT3);
    dist_to_player = u32(max(s32(dist_to_player) - s32(chunk_radius), 0));
    const u32 target_chunk_lod = get_chunk_lod_from_distance(dist_to_player);
    return target_chunk_lod;
}

struct GenChunkWork
{
    s32 chunk_x;
    s32 chunk_y;
    s32 chunk_z;
    u32 lod;
};
void allocate_and_gen_chunk_work(ThreadState* thread_state, void** out_data, const void* in_data)
{
    const GenChunkWork* data = reinterpret_cast<const GenChunkWork*>(in_data);
    Chunk* gen_chunk_scratch = thread_state->gen_chunk_scratch;
    assert(is_power_of_2(VOXELS_CHUNK_DIM));
    const s32 chunk_x = data->chunk_x;
    const s32 chunk_y = data->chunk_y;
    const s32 chunk_z = data->chunk_z;
    const s32 lod = data->lod;
    u32 num_voxels = 0;

    const s32 lod_scale = get_lod_scale(lod);
    assert(lod_scale != 0);
    
    // Generate a bitcube of which voxels are terrain.
    using ChunkBitCube = BitCube<VOXELS_CHUNK_DIM + 16, VOXELS_CHUNK_DIM + 2, VOXELS_CHUNK_DIM + 2>;
    ChunkBitCube bitcube;
    static_assert(sizeof(BitCube<VOXELS_CHUNK_DIM, VOXELS_CHUNK_DIM, VOXELS_CHUNK_DIM>) < 1024*1024);
    static_assert(sizeof(bitcube.m_v[0]) == 1);
    {
        const __m256 x_base = _mm256_add_ps(
                                            _mm256_set1_ps(f32(chunk_x - 8*lod_scale)),
                                            _mm256_mul_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(f32(lod_scale)))
                                            );
        __m256 z = _mm256_set1_ps(f32(chunk_z - f32(lod_scale)));
        for(s32 zi = 0; zi < VOXELS_CHUNK_DIM + 2; zi++)
        {
            __m256 y = _mm256_set1_ps(f32(chunk_y - lod_scale));
            for(s32 yi = 0; yi < VOXELS_CHUNK_DIM + 2; yi++)
            {
                __m256 x_accum = _mm256_set1_ps(0.0f);
                for(s32 xi = 0; xi < VOXELS_CHUNK_DIM + 16; xi += 8)
                {
                    __m256 x = _mm256_add_ps(x_base, x_accum);
                    __m256 n = sample_terrain(x, y, z);
                    __m256 n_mask = _mm256_cmp_ps(n, _mm256_set1_ps(0.0f), _CMP_NLT_UQ);
                    u8 bits = u8(_mm256_movemask_ps(n_mask));
                    u32 idx = (zi*(VOXELS_CHUNK_DIM+2)*(VOXELS_CHUNK_DIM+16)) + (yi*(VOXELS_CHUNK_DIM+16)) + xi;
                    bitcube.m_v[idx/8] = bits;
                    
                    x_accum = _mm256_add_ps(x_accum, _mm256_set1_ps(8.0f*f32(lod_scale)));
                }
                y = _mm256_add_ps(y, _mm256_set1_ps(f32(lod_scale)));
            }
            z = _mm256_add_ps(z, _mm256_set1_ps(f32(lod_scale)));
        }
    }
    
    ChunkBitCube bitcube_copy;
    memcpy(bitcube_copy.m_v, bitcube.m_v, ARRAY_COUNT(bitcube.m_v));
    
    // Turn off bits that are surrounded.
    for(s32 bitcube_zi = 1; bitcube_zi < VOXELS_CHUNK_DIM + 1; bitcube_zi++)
    {
        for(s32 bitcube_yi = 1; bitcube_yi < VOXELS_CHUNK_DIM + 1; bitcube_yi++)
        {
            for(s32 bitcube_xi = 8; bitcube_xi < VOXELS_CHUNK_DIM + 8; bitcube_xi++)
            {
                s32 idx = (bitcube_zi*(VOXELS_CHUNK_DIM+2)*(VOXELS_CHUNK_DIM+16)) + (bitcube_yi*(VOXELS_CHUNK_DIM+16)) + bitcube_xi;
                const u8 b0 = bitcube_copy.is_bit_set(bitcube_xi + 1, bitcube_yi + 0, bitcube_zi + 0);
                const u8 b1 = bitcube_copy.is_bit_set(bitcube_xi - 1, bitcube_yi + 0, bitcube_zi + 0);
                const u8 b2 = bitcube_copy.is_bit_set(bitcube_xi + 0, bitcube_yi + 1, bitcube_zi + 0);
                const u8 b3 = bitcube_copy.is_bit_set(bitcube_xi + 0, bitcube_yi - 1, bitcube_zi + 0);
                const u8 b4 = bitcube_copy.is_bit_set(bitcube_xi + 0, bitcube_yi + 0, bitcube_zi + 1);
                const u8 b5 = bitcube_copy.is_bit_set(bitcube_xi + 0, bitcube_yi + 0, bitcube_zi - 1);
                const u8 b6 = bitcube_copy.is_bit_set(bitcube_xi, bitcube_yi, bitcube_zi);
                
                const u8 mask = (b0<<0) | (b1<<1) | (b2<<2) | (b3<<3) | (b4<<4) | (b5<<5);
                
                bitcube.assign_bit(idx, (mask == 0b0011'1111) && b6 ? 0 : b6);
                                         
            }
        }
    }
    
    // TODO(mfritz) delete
    static u32 c = 4242;
    c = _mm_crc32_u32(c, chunk_z*MAX_CHUNKS*MAX_CHUNKS + chunk_y*MAX_CHUNKS + chunk_x);
    
    // Left pack and output.
    {
        const __m256 x_base = _mm256_add_ps(
                _mm256_set1_ps(f32(chunk_x - 8)),
                _mm256_mul_ps(_mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f), _mm256_set1_ps(f32(lod_scale)))
                );
        for(s32 bitcube_zi = 1; bitcube_zi < VOXELS_CHUNK_DIM+1; bitcube_zi++)
        {
            for(s32 bitcube_yi = 1; bitcube_yi < VOXELS_CHUNK_DIM+1; bitcube_yi++)
            {
                for(s32 bitcube_xi = 1; bitcube_xi < VOXELS_CHUNK_DIM/8 + 1; bitcube_xi++)
                {
                    u32 row_bytes = (VOXELS_CHUNK_DIM/8) + 2;
                    u32 frame_bytes = row_bytes*(VOXELS_CHUNK_DIM+2);
                    u32 bitcube_idx = bitcube_zi*frame_bytes + bitcube_yi*row_bytes + bitcube_xi;
                    u8 bits = bitcube.m_v[bitcube_idx];

                    u32 xi = (bitcube_xi - 1)*8;
                    u32 yi = (bitcube_yi - 1);
                    u32 zi = (bitcube_zi - 1);

                    __m256i x = _mm256_add_epi32(
                            _mm256_set1_epi32(chunk_x + xi*lod_scale),
                            _mm256_mullo_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), _mm256_set1_epi32(lod_scale))
                            );
                    __m256i y = _mm256_set1_epi32(chunk_y + yi*lod_scale);
                    __m256i z = _mm256_set1_epi32(chunk_z + zi*lod_scale);

                    u32 r = 0x22;
                    u32 g = 0x99;
                    u32 b = 0x00;
                    __m256 inten = pnoise8(
                            _mm256_mul_ps(_mm256_cvtepi32_ps(x), _mm256_set1_ps(0.01f)),
                            _mm256_mul_ps(_mm256_cvtepi32_ps(y), _mm256_set1_ps(0.01f)),
                            _mm256_mul_ps(_mm256_cvtepi32_ps(z), _mm256_set1_ps(0.01f))
                        );
                    __m256 vr = _mm256_add_ps(_mm256_set1_ps(f32(r)), _mm256_mul_ps(inten, _mm256_set1_ps(64.0f)));
                    __m256 vg = _mm256_add_ps(_mm256_set1_ps(f32(g)), _mm256_mul_ps(inten, _mm256_set1_ps(64.0f)));
                    __m256 vb = _mm256_add_ps(_mm256_set1_ps(f32(b)), _mm256_mul_ps(inten, _mm256_set1_ps(0.0f)));
                    __m256i vri = _mm256_cvtps_epi32(vr);
                    __m256i vgi = _mm256_cvtps_epi32(vg);
                    __m256i vbi = _mm256_cvtps_epi32(vb);
                    __m256i v = _mm256_or_si256(
                        _mm256_slli_epi32(vri, 24),
                        _mm256_or_si256(
                            _mm256_slli_epi32(vgi, 16),
                            _mm256_slli_epi32(vbi, 8)
                        )
                    );
                    

                    //__m256i v = _mm256_set1_epi32(0x22AA00FF);

                    /*
                    static constexpr u32 lod_lut[] = {
                        0x00AA00FF,
                        0x0000AAFF,
                        0xAA0000FF,
                        0x00AAAAFF,
                        0xAAAA00FF,
                        0xAAAAAAFF
                    };
                    __m256i v = _mm256_set1_epi32(lod_lut[min(lod, ARRAY_COUNT(lod_lut) - 1)]);
                    */

                    __m256i x_packed = left_pack(x, bits);
                    __m256i y_packed = left_pack(y, bits);
                    __m256i z_packed = left_pack(z, bits);
                    __m256i v_packed = left_pack(v, bits);
                    _mm256_storeu_si256((__m256i*)(gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*0 + num_voxels), x_packed);
                    _mm256_storeu_si256((__m256i*)(gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*1 + num_voxels), y_packed);
                    _mm256_storeu_si256((__m256i*)(gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*2 + num_voxels), z_packed);
                    _mm256_storeu_si256((__m256i*)(gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*3 + num_voxels), v_packed);
                    num_voxels += _mm_popcnt_u32(bits);
                }
            }
        }
    }
    
    Chunk* result = allocate_chunk(num_voxels);
    result->bl_pos = v3i(chunk_x, chunk_y, chunk_z);
    result->lod = lod;
    result->num = num_voxels;
    memcpy(result->voxels + num_voxels*0, gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*0, num_voxels * sizeof(result->voxels[0]));
    memcpy(result->voxels + num_voxels*1, gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*1, num_voxels * sizeof(result->voxels[0]));
    memcpy(result->voxels + num_voxels*2, gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*2, num_voxels * sizeof(result->voxels[0]));
    memcpy(result->voxels + num_voxels*3, gen_chunk_scratch->voxels + CHUNK_MAX_VOXELS*3, num_voxels * sizeof(result->voxels[0]));
    *out_data = result;
}



struct OctTree
{
    struct Node
    {
        Chunk* chunk;
        v3i bl_pos;
        u32 lod;
        u32 children[8];
    };
    static constexpr u32 CAP = MAX_CHUNKS;
    u32 num_allocated;
    Node pool[CAP];
    using NodeBitArray = BitArray<MAX_CHUNKS>;
    NodeBitArray leaves;

    Node* root()
    {
        return &(pool[0]);
    }

    const Node* root() const
    {
        return &(pool[0]);
    }

    u32 allocate(v3i pos, u32 lod)
    {
        u32 result_id = num_allocated++;
        assert(num_allocated < MAX_CHUNKS);
        Node* node = root() + result_id;
        node->bl_pos = pos;
        node->lod = lod;
        node->chunk = nullptr;
        return result_id;
    }

    u32 find(v3i bl_pos, u32 lod) const
    {
        const Node* cur = root();
        while(true)
        {
            const u32 cur_id = u32(cur - root());
            if(bl_pos == cur->bl_pos && lod == cur->lod)
            {
                return cur_id;
            }
            if(leaves.is_bit_set(cur_id) || lod > cur->lod)
            {
                return CAP;
            }

            bool progress = false;
            for(u32 i = 0; i < 8; i++)
            {
                // TODO Don't need to read from this memory if we statically know the order of quadrants...
                const Node* child = &(pool[cur->children[i]]);
                const v3i quadrant_bl = child->bl_pos;
                const u32 quadrant_dim = get_chunk_dim_from_lod(child->lod);
                if(point_cube_intersect(bl_pos, quadrant_bl, quadrant_dim))
                {
                    progress = true;
                    cur = child;
                    break;
                }
            }
            assert(progress);
            if(!progress)
            {
                return CAP;
            }
        }
    }

    bool is_leaf(u32 id) const
    {
        if(id == CAP) return false;
        return leaves.is_bit_set(id);
    }

    void reset()
    {
        // 1 node always allocated for the root.
        num_allocated = 1;
        leaves.clear();
        leaves.set_bit(0);
        const u32 root_node_hdim = get_chunk_dim_from_lod(MAX_LOD) / 2;
        root()->bl_pos = v3i(-s32(root_node_hdim), -s32(root_node_hdim), -s32(root_node_hdim));
        root()->lod = MAX_LOD;
    }
};
void dealloc_children(OctTree* ot, OctTree::Node* cur)
{
    for(u32 i = 0; i < 8; i++)
    {
        const u32 child_id = cur->children[i];
        OctTree::Node* child = ot->root() + child_id;
        deallocate_chunk(&child->chunk);
        if(!ot->is_leaf(child_id))
        {
            dealloc_children(ot, child);
        }
    }
};
void terrain_generation(
    OctTree* current_oct_tree,
    OctTree* target_oct_tree,
    const v3i player_pos)
{
    // Build the oct tree for terrain.
    u32* nodes_to_process = new u32[MAX_CHUNKS];
    u32 num_nodes_to_process = 0;
    {
        TIME_SCOPE("Build oct tree");

        nodes_to_process[num_nodes_to_process++] = 0;

        // Process until all nodes are the appropriate size.
        while(num_nodes_to_process)
        {
            // Pop next node.
            u32 node_id = nodes_to_process[--num_nodes_to_process];
            assert(node_id < OctTree::CAP);
            OctTree::Node* node = target_oct_tree->pool + node_id;

            // Do split?
            const u32 target_chunk_lod = get_chunk_target_lod(node->bl_pos, node->lod, player_pos);
            if(node->lod > target_chunk_lod)
            {
                // Split
                const u32 new_lod = node->lod - 1;
                assert(new_lod < MAX_LOD);
                const u32 new_dim = get_chunk_dim_from_lod(new_lod);
                target_oct_tree->leaves.clear_bit(node_id);
                v3i child_offsets[8] = 
                {
                    v3i(0, 0, 0),
                    v3i(1, 0, 0),
                    v3i(0, 1, 0),
                    v3i(1, 1, 0),
                    v3i(0, 0, 1),
                    v3i(1, 0, 1),
                    v3i(0, 1, 1),
                    v3i(1, 1, 1),
                };
                for(u32 i = 0; i < 8; i++)
                {
                    u32 child_id = target_oct_tree->allocate(node->bl_pos + child_offsets[i] * new_dim, new_lod);
                    node->children[i] = child_id;
                    target_oct_tree->leaves.set_bit(child_id);

                    nodes_to_process[num_nodes_to_process++] = child_id;
                    assert(num_nodes_to_process < MAX_CHUNKS);
                }
            }
        }
    }


    u32 num_work = 0;
    GenChunkWork* work = new GenChunkWork[MAX_CHUNKS_GEN_PER_FRAME];
    OctTree::Node** work_output_nodes = new OctTree::Node*[MAX_CHUNKS_GEN_PER_FRAME];


    num_nodes_to_process = 1;
    nodes_to_process[0] = 0;
    while(num_nodes_to_process)
    {
        // Pop next node.
        u32 current_node_id = nodes_to_process[--num_nodes_to_process];
        assert(current_node_id < OctTree::CAP);
        OctTree::Node* current_node = current_oct_tree->root() + current_node_id;

        // Find current_node in new oct tree.
        u32 target_node_id = target_oct_tree->find(current_node->bl_pos, current_node->lod);
        assert(target_node_id != OctTree::CAP);

        // Combine
        if(!current_oct_tree->is_leaf(current_node_id) && target_oct_tree->is_leaf(target_node_id))
        {
            // If there is room to generate 1 chunk:
            // Generate the new parent chunk. All children will need to be deallocated.
            
            if(num_work < MAX_CHUNKS_GEN_PER_FRAME)
            {
                dealloc_children(current_oct_tree, current_oct_tree->root() + current_node_id);

                work[num_work] = GenChunkWork{
                    .chunk_x = current_node->bl_pos.x(),
                    .chunk_y = current_node->bl_pos.y(),
                    .chunk_z = current_node->bl_pos.z(),
                    .lod = current_node->lod
                };
                work_output_nodes[num_work] = current_oct_tree->root() + current_node_id;
                num_work++;

                current_oct_tree->leaves.set_bit(current_node_id);
            }
        }
        // Split
        else if(current_oct_tree->is_leaf(current_node_id) && !target_oct_tree->is_leaf(target_node_id))
        {
            // If there is room to generate 8 chunks:
            // Generate the 8 child chunks. Deallocate old parent chunk.

            if(num_work + 8 <= MAX_CHUNKS_GEN_PER_FRAME)
            {
                deallocate_chunk(&current_node->chunk);

                // Only allocate immediate children for now...
                for(u32 i = 0; i < 8; i++)
                {
                    OctTree::Node* target_node = target_oct_tree->root() + target_node_id;
                    OctTree::Node* target_child = target_oct_tree->root() + target_node->children[i];
                    work[num_work] = GenChunkWork{
                        .chunk_x = target_child->bl_pos.x(),
                        .chunk_y = target_child->bl_pos.y(),
                        .chunk_z = target_child->bl_pos.z(),
                        .lod = target_child->lod
                    };

                    u32 current_child_id = current_oct_tree->allocate(target_child->bl_pos, target_child->lod);
                    current_node->children[i] = current_child_id;
                    current_oct_tree->leaves.set_bit(current_child_id);

                    work_output_nodes[num_work] = current_oct_tree->root() + current_child_id;
                    num_work++;
                }

                current_oct_tree->leaves.clear_bit(current_node_id);
            }
        }
        // Continue traversal
        else if(!current_oct_tree->is_leaf(current_node_id) && !target_oct_tree->is_leaf(target_node_id))
        {
            for(u32 i = 0; i < 8; i++)
            {
                nodes_to_process[num_nodes_to_process++] = current_node->children[i];
            }
        }
    }

    for(u32 i = 0; i < num_work; i++)
    {
        add_work(reinterpret_cast<void**>(&work_output_nodes[i]->chunk), JobId::load_terrain, allocate_and_gen_chunk_work, work + i);
    }

    wait_for_job(JobId::load_terrain);

    delete[] nodes_to_process;
    nodes_to_process = nullptr;

    delete[] work_output_nodes;
    work_output_nodes = nullptr;

    delete[] work;
    work = nullptr;
}


void voxel_line(
        s32* out_pos,
        u32* out_voxel_id,
        u32* out_num,
        const u32 pos_stride,
        const v3 a,
        const v3 b,
        const f32 a_thickness,
        const f32 b_thickness)
{
    f32 max_thickness = max(a_thickness, b_thickness);
    s32 min_x = min(s32(a.x), s32(b.x)) - s32(ceil(max_thickness));
    s32 min_y = min(s32(a.y), s32(b.y)) - s32(ceil(max_thickness));
    s32 min_z = min(s32(a.z), s32(b.z)) - s32(ceil(max_thickness));
    s32 max_x = max(s32(a.x), s32(b.x)) + s32(ceil(max_thickness));
    s32 max_y = max(s32(a.y), s32(b.y)) + s32(ceil(max_thickness));
    s32 max_z = max(s32(a.z), s32(b.z)) + s32(ceil(max_thickness));
    
    v3 d = b - a;
    f32 d_len = sqrtf(dot(d, d));
    v3 dn = d / d_len;
    u32 num = *out_num;
    for(s32 z = min_z; z <= max_z; z++)
    {
        for(s32 y = min_y; y <= max_y; y++)
        {
            for(s32 x = min_x; x <= max_x; x++)
            {
                v3 p = v3(f32(x), f32(y), f32(z));
                v3 v = p - a;
                f32 t = dot(v, dn);
                v3 n = (p - dn*t) - a;
                f32 dist = sqrtf(n.x*n.x + n.y*n.y + n.z*n.z);
                f32 tn = t / d_len;
                f32 thickness = lerp(a_thickness, b_thickness, tn);
                dist = t < 0.0f ? length(a - p) : dist;
                dist = t > d_len ? length(b - p) : dist;
                if(dist < thickness)
                {
                    out_pos[pos_stride*0 + num] = x;
                    out_pos[pos_stride*1 + num] = y;
                    out_pos[pos_stride*2 + num] = z;
                    out_voxel_id[num] = 0x442222FF;
                    num++;
                }
            }
        }
    }
    *out_num = num;
};

struct TreeParams
{
    v3 base = v3();
    v3 dir = v3();
    f32 len = 0.0f;
    u32 count = 0;
    f32 thickness = 0.0f;
    f32 thickness_falloff = 0.0f;
    f32 variance = 0.0f;
};
void make_tree_branch(Tree& out, const TreeParams& in)
{
    ASSERT(in.count < Tree::MAX_NODES, "Too many nodes for tree. Must be less than Tree::MAX_NODES - 1");
    
    f32 len = in.len;
    f32 variance = in.variance;
    v3 dir = in.dir;
    
    u32 num_nodes = 0;
    v3 cur_node_pos = in.base;
    
    out.nodes_pos[num_nodes] = cur_node_pos;
    out.nodes_dist[num_nodes] = 0;
    num_nodes++;
    
    f32 len_accum = 0.0f;
    for(u32 i = 0; i < in.count; i++)
    {
        cur_node_pos += dir * len;
        
        out.nodes_pos[num_nodes] = cur_node_pos;
        out.nodes_dist[num_nodes] = len_accum;
        num_nodes++;
        
        len_accum += len;
        
        v3 i_axis = dir;
        v3 j_axis = cross(dir, v3(0.0f, 1.0f, 0.0f));
        v3 k_axis = cross(i_axis, j_axis);
        dir += random_range(-variance, variance) * j_axis + random_range(-variance, variance) * k_axis;
        dir += v3(0.0f, variance * 0.1f, 0.0f); // Push tree upwards
        dir = normalize(dir);
        len *= 0.9f;
        variance *= 1.1f;
    }
    
    out.num_nodes = num_nodes;
}
void make_tree(Tree& out, v3 base)
{
    make_tree_branch(out, TreeParams{
                         .base = base,
                         .dir = normalize(v3(0.0f, 1.0f, 0.1f)),
                         .len = 12.0f,
                         .count = 6,
                         .thickness = 0.0f,
                         .thickness_falloff = 0.0f,
                         .variance = 2.0f,
                     });
    f32 h = 12.0f;
    out.num_branches = 0;
    for(u32 i = 0; i < 4; i++)
    {
        // TODO(mfritz) No dynamic allocation
        out.branches[i] = new Tree;
        make_tree_branch(*out.branches[i], TreeParams{
                             .base = out.get_point_at_height(h),
                             .dir = make_axis_rotation_matrix(random_range(0.0f, 2*M_PI), v3(0.0f, 1.0f, 0.0f)) * normalize(v3(0.0f, 1.0f, 1.0f)),
                             .len = 5.0f,
                             .count = 4,
                             .thickness = 0.0f,
                             .thickness_falloff = 0.0f,
                             .variance = 0.5f,
                         });
        out.num_branches++;
        
        h += 10.0f;
    }
}
void leaf_ball(s32* out_pos, u32* out_voxel_id, u32* out_num, const u32 pos_stride, const v3 leaf_center, s32 radius)
{
    u32 num = *out_num;
    ITER_CUBE(radius*2)
    {
        s32 px = it.x - radius + s32(floor(leaf_center.x));
        s32 py = it.y - radius + s32(floor(leaf_center.y));
        s32 pz = it.z - radius + s32(floor(leaf_center.z));
        v3 p = v3(float(px), float(py), float(pz));
        //f32 r = radius - perlin_01(p.x * 20.0f, p.y * 20.0f, p.z * 20.0f) * (f32(radius) * 0.5f);
        f32 r = radius - 0.5f * f32(radius) * 0.5f;
        if(abs(length(p - leaf_center) - r) < 1.5f)
        {
            out_pos[pos_stride*0 + num] = px;
            out_pos[pos_stride*1 + num] = py;
            out_pos[pos_stride*2 + num] = pz;
            out_voxel_id[num] = 0x50BA2CFF;
            num++;
        }
    }
    *out_num = num;
}

void rasterize_tree(s32* out_pos, u32* out_voxel_id, u32* out_num, const u32 pos_stride, const Tree& tree)
{
    for(u32 i = 0; i < tree.num_nodes - 1; i++)
    {
        voxel_line(out_pos, out_voxel_id, out_num, pos_stride, tree.nodes_pos[i], tree.nodes_pos[i+1], 2.0f, 2.0f);
        ASSERT(*out_num < TreeBuffer::CAP, "Rasterize tree too big");
    }
    
    leaf_ball(out_pos, out_voxel_id, out_num, pos_stride, tree.nodes_pos[tree.num_nodes - 1], 22);
    ASSERT(*out_num < TreeBuffer::CAP, "Rasterize tree too big");
    
    for(u32 b = 0; b < tree.num_branches; b++)
    {
        Tree *branch = tree.branches[b];
        for(u32 i = 0; i < branch->num_nodes - 1; i++)
        {
            voxel_line(out_pos, out_voxel_id, out_num, pos_stride, branch->nodes_pos[i], branch->nodes_pos[i+1], 1.0f, 1.0f);
            ASSERT(*out_num < TreeBuffer::CAP, "Rasterize tree too big");
        }
        leaf_ball(out_pos, out_voxel_id, out_num, pos_stride, branch->nodes_pos[branch->num_nodes - 1], 14);
        ASSERT(*out_num < TreeBuffer::CAP, "Rasterize tree too big");
    }
}

void debug_draw_tree(const Tree& tree)
{
    for(u32 i = 0; i < tree.num_nodes - 1; i++)
    {
        draw_line(tree.nodes_pos[i], tree.nodes_pos[i+1], v3(0.0f, 1.0f, 0.0f));
    }
    for(u32 b = 0; b < tree.num_branches; b++)
    {
        Tree *branch = tree.branches[b];
        for(u32 i = 0; i < branch->num_nodes - 1; i++)
        {
            draw_line(branch->nodes_pos[i], branch->nodes_pos[i+1], v3(0.0f, 1.0f, 0.0f));
        }
    }
}


INT WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR lpCmdLine, INT nCmdShow)
{
    (void)nCmdShow;
    (void)lpCmdLine;
    (void)hPrevInstance;
    (void)hInstance;

    G_log_file = fopen("log.txt", "wt");
    
    main_thread_state = new ThreadState();
    main_thread_state->gen_chunk_scratch = allocate_chunk(CHUNK_MAX_VOXELS);
    if constexpr (NUM_TOTAL_THREADS > 1)
    {
        constexpr u32 NUM_NEW_THREADS = NUM_TOTAL_THREADS - 1;
        // + 1 so the compiler won't complain when total thread count is 1
        HANDLE threads[NUM_NEW_THREADS + 1];
        u32 thread_ids[NUM_NEW_THREADS + 1];
        for(u32 i = 0; i < NUM_NEW_THREADS; i++)
        {
            ThreadState* thread_state = new ThreadState();
            thread_state->gen_chunk_scratch = allocate_chunk(CHUNK_MAX_VOXELS);
            threads[i] = CreateThread( 
                NULL,                   // default security attributes
                0,                      // use default stack size  
                worker_main,            // thread function name
                thread_state,           // argument to thread function 
                0,                      // use default creation flags 
                reinterpret_cast<DWORD*>(&thread_ids[i])); // returns the thread identifier 
            if(threads[i] == NULL)
            {
                return 1;
            }
            
            u32 affinity_mask = 1 << i;
            DWORD_PTR old_affinity_mask = SetThreadAffinityMask(threads[i], affinity_mask);
            if(old_affinity_mask == 0)
            {
                return 1;
            }
            SetThreadDescription(threads[i], L"worker_thread");
        }
    }
    
    // GLFW
    G_graphics_state = new GraphicsState();
    {
        glfwSetErrorCallback(glfw_error_fn);
        if(!glfwInit())
        {
            fprintf(G_log_file, "Failed to initialize GLFW\n");
            return 1;
        }
        
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_SAMPLES, 4); // TODO Not working?
        
        //u32 window_width = 1920 * 1.0f;
        //u32 window_height = 1080 * 1.0f;
        u32 window_width = 2560;
        u32 window_height = 1440;
        G_graphics_state->window = glfwCreateWindow(window_width, window_height, "My Game", NULL, NULL);
        if(!G_graphics_state->window)
        {
            fprintf(G_log_file, "Failed to open GLFW window\n");
            glfwTerminate();
            return 1;
        }
        
        glfwSetFramebufferSizeCallback(G_graphics_state->window, reshape);
        //glfwSetKeyCallback(G_graphics_state->window, key);
        
        glfwMakeContextCurrent(G_graphics_state->window);
        
        if(!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
        {
            fprintf(G_log_file, "Failed to initialize OpenGL context\n");
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
        glfwGetFramebufferSize(G_graphics_state->window, &framebuffer_width, &framebuffer_height);
        reshape(G_graphics_state->window, framebuffer_width, framebuffer_height);
        
        G_graphics_state->cam.pos = v3(16.0f, 16.0f, 16.0f);
        G_graphics_state->cam.rot_x = 0.0f;
        G_graphics_state->cam.rot_y = 0.0f;
        G_graphics_state->cam.vfov = 90.0f;
        G_graphics_state->cam.view_dist = 5000.0f;
        G_graphics_state->cam.near_plane_dist = 1.0f;
        G_graphics_state->debug_cam = G_graphics_state->cam;
        
        {
            glGenVertexArrays(1, &G_graphics_state->batch_voxel_vao);
            glBindVertexArray(G_graphics_state->batch_voxel_vao);
            check_gl_errors("vao");
            glGenBuffers(1, &G_graphics_state->batch_voxel_vbo);
            check_gl_errors("vbo");
            
            glBindBuffer(GL_ARRAY_BUFFER, G_graphics_state->batch_voxel_vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(GraphicsState::voxel_vertices), GraphicsState::voxel_vertices, GL_STATIC_DRAW);
            check_gl_errors("vbo voxel vertices");
            
            glGenBuffers(1, &G_graphics_state->batch_voxel_ebo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, G_graphics_state->batch_voxel_ebo);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GraphicsState::voxel_indices), GraphicsState::voxel_indices, GL_STATIC_DRAW);
            
            glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void *)0);
            glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (void *)(8*sizeof(f32)));
            glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, (void *)(16*sizeof(f32)));
            glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, 0, (void *)(24*sizeof(f32)));
            glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, 0, (void *)(32*sizeof(f32)));
            glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, 0, (void *)(40*sizeof(f32)));
            glEnableVertexAttribArray(0);
            glEnableVertexAttribArray(1);
            glEnableVertexAttribArray(2);
            glEnableVertexAttribArray(3);
            glEnableVertexAttribArray(4);
            glEnableVertexAttribArray(5);
            check_gl_errors("vertex attrib pointer");
            
            // https://www.khronos.org/opengl/wiki/Shader_Storage_Buffer_Object
            const u32 voxel_buffer_size = SHADER_BUFFER_WIDTH * Chunk::INTS_PER_VOXEL * 5;
            glGenBuffers(1, &G_graphics_state->batch_voxel_ssbo);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, G_graphics_state->batch_voxel_ssbo);
            glBufferData(GL_SHADER_STORAGE_BUFFER, voxel_buffer_size, nullptr, GL_DYNAMIC_DRAW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, G_graphics_state->batch_voxel_ssbo);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            check_gl_errors("ssbo");
            
            G_graphics_state->batch_voxel_shader_program = make_shader_from_file("assets/shaders/batch_voxel.shader");
            {
                fprintf(G_log_file, G_graphics_state->debug_info_log);
            }
        }
        
        {
            glGenVertexArrays(1, &G_graphics_state->debug_line_vao);
            glBindVertexArray(G_graphics_state->debug_line_vao);
            check_gl_errors("vao");
            glGenBuffers(1, &G_graphics_state->debug_line_vbo);
            check_gl_errors("vbo");
            glBindBuffer(GL_ARRAY_BUFFER, G_graphics_state->debug_line_vbo);
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
            
            G_graphics_state->debug_line_shader_program = make_shader_from_string(
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
                fprintf(G_log_file, G_graphics_state->debug_info_log);
                fclose(file);
            }
        }
        
        glGenTextures(1, &G_graphics_state->imgui_debug_texture);
        glBindTexture(GL_TEXTURE_2D, G_graphics_state->imgui_debug_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        
        G_graphics_state->imgui_debug_texture_width = 512;
        G_graphics_state->imgui_debug_texture_height = 512;
        G_graphics_state->imgui_debug_texture_data = new u32[G_graphics_state->imgui_debug_texture_width * G_graphics_state->imgui_debug_texture_height]{};
        glTexImage2D(
                     GL_TEXTURE_2D,
                     0,
                     GL_RGBA,
                     G_graphics_state->imgui_debug_texture_width,
                     G_graphics_state->imgui_debug_texture_height,
                     0,
                     GL_RGBA,
                     GL_UNSIGNED_BYTE,
                     G_graphics_state->imgui_debug_texture_data);
    }
    
    // ImGui
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplOpenGL3_Init("#version 440 core");
        ImGui_ImplGlfw_InitForOpenGL(G_graphics_state->window, true);
    }
    
    G_input_state = new InputState();
    
    PackedVoxels* packed_voxels = new PackedVoxels;
    VoxelRenderData * voxel_render_data = new VoxelRenderData;

    //u32 num_chunks = 0;
    //Chunk** chunks = new Chunk*[MAX_CHUNKS];

    OctTree* oct_tree = new OctTree;
    oct_tree->reset();

    OctTree* target_oct_tree = new OctTree;
    target_oct_tree->reset();
    
#if 0
    static constexpr u32 NUM_TREES = 1;
    Tree trees[NUM_TREES];
    const u32 num_trees = NUM_TREES;
    f32 tree_pos[NUM_TREES*3];
    for(u32 i = 0; i < NUM_TREES; i++)
    {
        make_tree(trees[i], v3());
        tree_pos[NUM_TREES*0 + i] = random_range(-128.0f, 128.0f);
        //tree_pos[NUM_TREES*1 + i] = random_range(-32.0f, 32.0f);
        tree_pos[NUM_TREES*1 + i] = random_range(0.0f, 2.0f);
        tree_pos[NUM_TREES*2 + i] = random_range(-128.0f, 128.0f);
    }
    
    TreeBuffer* tree_buffers[NUM_TREES];
    for(u32 i = 0; i < ARRAY_COUNT(tree_buffers); i++)
    {
        // TODO(mfritz) No dynamic allocation.
        tree_buffers[i] = new TreeBuffer;
        tree_buffers[i]->num = 0;
        rasterize_tree(tree_buffers[i]->pos, tree_buffers[i]->voxel_id, &tree_buffers[i]->num, TreeBuffer::CAP, trees[i]);
    }
#endif
    
    f32 frame_timer = 0.0f;
    f32 last_time = 0.0f;
    const static f32 TIME_STEP = 0.016f;
    b32 running = true;
    while(running)
    {
        f64 current_time = glfwGetTime();
        frame_timer += f32(current_time) - f32(last_time);
        last_time = f32(current_time);
        
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
            
            const ImGuiTabBarFlags imgui_tab_bar_flags = ImGuiTabBarFlags_None;
            const bool imgui_tab_bar = ImGui::BeginTabBar("MyTabBar", imgui_tab_bar_flags);
            
            running = !glfwGetKey(G_graphics_state->window, GLFW_KEY_ESCAPE) && !glfwWindowShouldClose(G_graphics_state->window);
            
            static f32 camera_speed = 100.0f;
            if(imgui_tab_bar && ImGui::BeginTabItem("Debug Camera"))
            {
                ImGui::InputFloat("camera x", &G_graphics_state->cam.pos.x, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::InputFloat("camera y", &G_graphics_state->cam.pos.y, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::InputFloat("camera z", &G_graphics_state->cam.pos.z, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::Text("camera x rot %f", G_graphics_state->cam.rot_x);
                ImGui::Text("camera y rot %f", G_graphics_state->cam.rot_y);
                ImGui::InputFloat("debug camera x", &G_graphics_state->debug_cam.pos.x, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::InputFloat("debug camera y", &G_graphics_state->debug_cam.pos.y, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::InputFloat("debug camera z", &G_graphics_state->debug_cam.pos.z, 1.0f, 10.0f, 3, ImGuiInputTextFlags_EnterReturnsTrue);
                ImGui::DragFloat("camera near_plane_dist", &G_graphics_state->cam.near_plane_dist);
                ImGui::DragFloat("camera view_dist", &G_graphics_state->cam.view_dist);
                ImGui::Checkbox("debug cam", &G_graphics_state->use_debug_cam);
                
                static float d_move[3] = {};
                ImGui::InputFloat("move by x", &d_move[0], 1.0f, 10.0f, 3);
                ImGui::InputFloat("move by y", &d_move[1], 1.0f, 10.0f, 3);
                ImGui::InputFloat("move by z", &d_move[2], 1.0f, 10.0f, 3);
                if(ImGui::Button("move"))
                {
                    G_graphics_state->cam.pos.x += d_move[0];
                    G_graphics_state->cam.pos.y += d_move[1];
                    G_graphics_state->cam.pos.z += d_move[2];
                }
                ImGui::DragFloat("camera speed", &camera_speed);
                ImGui::EndTabItem();
            }
            
            {
                TIME_SCOPE("move camera");
                
                Camera *current_cam;
                if(G_graphics_state->use_debug_cam) current_cam = &G_graphics_state->debug_cam;
                else current_cam = &G_graphics_state->cam;
                
                static f32 last_mouse_x, last_mouse_y;
                f64 xpos, ypos;
                glfwGetCursorPos(G_graphics_state->window, &xpos, &ypos);
                if(glfwGetKey(G_graphics_state->window, GLFW_KEY_SPACE))
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
                if(glfwGetKey(G_graphics_state->window, GLFW_KEY_S)) camera_vel += -cam_i;
                if(glfwGetKey(G_graphics_state->window, GLFW_KEY_F)) camera_vel +=  cam_i;
                if(glfwGetKey(G_graphics_state->window, GLFW_KEY_D)) camera_vel +=  cam_k;
                if(glfwGetKey(G_graphics_state->window, GLFW_KEY_E)) camera_vel += -cam_k;
                current_cam->pos += normalize(camera_vel) * TIME_STEP * camera_speed;
                
            }
            
            const v3i camera_pos = point_to_voxel(v3(G_graphics_state->cam.pos.x, G_graphics_state->cam.pos.y, G_graphics_state->cam.pos.z));

            // Gen terrain
            {
                TIME_SCOPE("Generate terrain");

                target_oct_tree->reset();
                terrain_generation(oct_tree, target_oct_tree, camera_pos);
            }

            // Prep render data
#if 1
            {
                TIME_SCOPE("Prep render data");

                u32 num_voxels = 0;
                u32* nodes_to_process = new u32[MAX_CHUNKS];
                u32 num_nodes_to_process = 1;
                nodes_to_process[0] = 0;
                while(num_nodes_to_process)
                {
                    // Pop next node.
                    u32 node_id = nodes_to_process[--num_nodes_to_process];
                    OctTree::Node* node = oct_tree->root() + node_id;

                    if(node->chunk)
                    {
                        const Chunk* chunk = node->chunk;
                        memcpy(packed_voxels->pos + MAX_VOXELS*0 + num_voxels, chunk->voxels + chunk->num*0, chunk->num * sizeof(chunk->voxels[0]));
                        memcpy(packed_voxels->pos + MAX_VOXELS*1 + num_voxels, chunk->voxels + chunk->num*1, chunk->num * sizeof(chunk->voxels[0]));
                        memcpy(packed_voxels->pos + MAX_VOXELS*2 + num_voxels, chunk->voxels + chunk->num*2, chunk->num * sizeof(chunk->voxels[0]));
                        memcpy(packed_voxels->voxel_id + num_voxels, chunk->voxels + chunk->num*3, chunk->num * sizeof(chunk->voxels[0]));
                        
                        for(u32 i = num_voxels; i < num_voxels + chunk->num; i++)
                        {
                            packed_voxels->scale[i] = 1 << chunk->lod;
                        }
                        
                        num_voxels += chunk->num;
                    }

                    if(!oct_tree->is_leaf(node_id))
                    {
                        for(u32 i = 0; i < 8; i++)
                        {
                            assert(node->children[i] < OctTree::CAP);
                            nodes_to_process[num_nodes_to_process++] = node->children[i];
                        }
                    }
                }
                packed_voxels->num = num_voxels;
                
                delete[] nodes_to_process;
                nodes_to_process = nullptr;
                
                // TODO(mfritz) Pad to nearest next 8 voxels with 0xFFFFFFFF
                // TODO(mfritz) Do LOD
                camera_cull(
                            &voxel_render_data->num,
                            voxel_render_data->pos,
                            voxel_render_data->scale,
                            voxel_render_data->color,
                            packed_voxels->num,
                            packed_voxels->pos,
                            packed_voxels->scale,
                            packed_voxels->voxel_id);
            }
#endif
            
            
            // Draw
            u32 num_batches_drawn = 0;
            {
                TIME_SCOPE("draw");
                
                for(u32 batch_i = 0; batch_i < voxel_render_data->num; batch_i += VoxelRenderData::BATCH_SIZE)
                {
                    f32 *batch_x = voxel_render_data->pos + MAX_VOXELS*0 + batch_i;
                    f32 *batch_y = voxel_render_data->pos + MAX_VOXELS*1 + batch_i;
                    f32 *batch_z = voxel_render_data->pos + MAX_VOXELS*2 + batch_i;
                    f32 *batch_scale = voxel_render_data->scale + batch_i;
                    u32 *batch_color = voxel_render_data->color + batch_i;
                    
                    u32 batch_size = min(voxel_render_data->num - batch_i, VoxelRenderData::BATCH_SIZE);
                    
                    assert(VoxelRenderData::BATCH_SIZE <= SHADER_BUFFER_WIDTH);
                    glUseProgram(G_graphics_state->batch_voxel_shader_program);
                    check_gl_errors("use program");
                    
                    Camera *current_cam;
                    if(G_graphics_state->use_debug_cam) current_cam = &G_graphics_state->debug_cam;
                    else current_cam = &G_graphics_state->cam;
                    
                    {
                        mat4 m = view_m_world(current_cam);
                        GLint loc = glGetUniformLocation(G_graphics_state->batch_voxel_shader_program, "m_view");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }
                    {
                        mat4 m = clip_m_view(current_cam);
                        GLint loc = glGetUniformLocation(G_graphics_state->batch_voxel_shader_program, "m_proj");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }
                    
                    glBindVertexArray(G_graphics_state->batch_voxel_vao);
                    glBindBuffer(GL_ARRAY_BUFFER, G_graphics_state->batch_voxel_vbo);
                    u32 offset_x = 0 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_y = 1 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_z = 2 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_scale = 3 * SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_color = 4 * SHADER_BUFFER_WIDTH*sizeof(u32);
                    u32 size = sizeof(f32) * batch_size;
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, G_graphics_state->batch_voxel_ssbo);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_x, size, batch_x);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_y, size, batch_y);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_z, size, batch_z);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_scale, size, batch_scale);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_color, size, batch_color);
                    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, G_graphics_state->batch_voxel_ssbo);
                    
                    const u32 num_indices = ARRAY_COUNT(GraphicsState::voxel_indices);
                    glDrawElementsInstanced(GL_TRIANGLES, num_indices, GL_UNSIGNED_INT, 0, batch_size);
                    check_gl_errors("draw");
                    
                    num_batches_drawn++;
                }
            }
            
            
            // Stats
#if 1
            if(imgui_tab_bar && ImGui::BeginTabItem("Stats"))
            {
                //ImGui::Text("View dist: %i", VIEW_DIST);
                ImGui::Text("Max chunks: %i", MAX_CHUNKS);
                //ImGui::Text("Chunks table mem: %iKB", KB(MAX_CHUNKS * sizeof(Chunk*)));
                //u32 num_generated_chunks = 0;
                //u32 num_populated_chunks = 0;
                //u32 avg_voxels_per_populated_chunk = 0;
                //u32 max_voxels_per_chunk = 0;
                //u32 voxels_bytes = 0;
                //for(u32 i = 0; i < num_chunks; i++)
                //{
                //    if(chunks[i])
                //    {
                //        num_generated_chunks++;
                //        num_populated_chunks += chunks[i]->num > 0;
                //        avg_voxels_per_populated_chunk += chunks[i]->num;
                //        max_voxels_per_chunk = max(max_voxels_per_chunk, chunks[i]->num);
                //        voxels_bytes += sizeof(chunks[i]->num) + chunks[i]->num * sizeof(chunks[i]->voxels[0]);
                //    }
                //}
                //ImGui::Text("%i / %i (%.1f%%) chunks generated", num_generated_chunks, MAX_CHUNKS, (f32(num_generated_chunks) / f32(MAX_CHUNKS)) * 100.0f);
                //ImGui::Text("# populated chunks: %i", num_populated_chunks);
                //ImGui::Text("# avg voxels per populated chunk: %i", num_populated_chunks == 0 ? 0 : avg_voxels_per_populated_chunk / num_populated_chunks);
                //ImGui::Text("# max voxels per chunk: %i", max_voxels_per_chunk);
                //ImGui::Text("Voxels mem: %i KB", KB(voxels_bytes));
                
                ImGui::Text("num_voxels %i", voxel_render_data->num);
                ImGui::Text("batches %i", num_batches_drawn);
                ImGui::Text("prep num_voxels %i", packed_voxels->num);
                
                ImGui::EndTabItem();
            }
#endif
            
            
            if(imgui_tab_bar && ImGui::BeginTabItem("Debug Draw"))
            {
                debug_draw_frustum();

                static bool show_chunk_lines = false;
                ImGui::Checkbox("show chunk links", &show_chunk_lines);
                static u32 show_lod = 0;
                ImGui::InputInt("show lod", reinterpret_cast<int*>(&show_lod));
                if(show_chunk_lines)
                {
                    // NOTE: Assuming oct tree nodes are allocated linearly and not deleted.
                    for(u32 node_id = 0; node_id < oct_tree->num_allocated; node_id++)
                    {
                        v3i bl_pos = oct_tree->pool[node_id].bl_pos;
                        u32 lod = oct_tree->pool[node_id].lod;

                        if(lod != show_lod) continue;

                        s32 dim = s32(get_chunk_dim_from_lod(lod));
                        v3i vs[] = 
                        {
                            bl_pos + v3i(  0,   0,   0),
                            bl_pos + v3i(dim,   0,   0),
                            bl_pos + v3i(  0, dim,   0),
                            bl_pos + v3i(dim, dim,   0),
                            bl_pos + v3i(  0,   0, dim),
                            bl_pos + v3i(dim,   0, dim),
                            bl_pos + v3i(  0, dim, dim),
                            bl_pos + v3i(dim, dim, dim) 
                        };
                        v3 vs_f[8];
                        for(u32 i = 0; i < 8; i++)
                        {
                            vs_f[i] = v3(f32(vs[i].x()), f32(vs[i].y()), f32(vs[i].z()));
                        }
                        f32 inten = remap(0.0f, 3.0f, 1.0f, 0.0f, f32(lod));
                        draw_line(vs_f[0], vs_f[1], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[1], vs_f[3], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[3], vs_f[2], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[2], vs_f[0], v3(1.0f, 0.0f, 1.0f) * inten);

                        draw_line(vs_f[4], vs_f[5], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[5], vs_f[7], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[7], vs_f[6], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[6], vs_f[4], v3(1.0f, 0.0f, 1.0f) * inten);

                        draw_line(vs_f[0], vs_f[4], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[1], vs_f[5], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[2], vs_f[6], v3(1.0f, 0.0f, 1.0f) * inten);
                        draw_line(vs_f[3], vs_f[7], v3(1.0f, 0.0f, 1.0f) * inten);
                    }
                }
                
                ImGui::EndTabItem();
            }
            
            
            // ImGui debug texture
            if(imgui_tab_bar && ImGui::BeginTabItem("Debug Texture"))
            {
                TIME_SCOPE("debug draw texture");
                memset(G_graphics_state->imgui_debug_texture_data, 0, G_graphics_state->imgui_debug_texture_width*G_graphics_state->imgui_debug_texture_height*sizeof(u32));
                ASSERT((G_graphics_state->imgui_debug_texture_width*G_graphics_state->imgui_debug_texture_height & 0b111) == 0, "Texture must be divisible by 8");
                
                static f32 x_off = 0.0f;
                static f32 y_off = 0.0f;
                static f32 z_off = 0.0f;
                ImGui::DragFloat("x_off", &x_off);
                ImGui::DragFloat("y_off", &y_off);
                ImGui::DragFloat("z_off", &z_off);
                ASSERT((G_graphics_state->imgui_debug_texture_width & 0b111) == 0, "Texture width must be divisible by 8");
                
                f32 texture_width = f32(G_graphics_state->imgui_debug_texture_width);
                f32 texture_height = f32(G_graphics_state->imgui_debug_texture_height);
                const __m256 x_base = _mm256_add_ps(_mm256_set1_ps(x_off - texture_width/2.0f), _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f));
                
                __m256 z = _mm256_set1_ps(z_off);
                __m256 y = _mm256_set1_ps(y_off - texture_height/2.0f);
                for(u32 yi = 0; yi < G_graphics_state->imgui_debug_texture_height; yi++)
                {
                    __m256 x_accum = _mm256_set1_ps(0.0f);
                    for(u32 xi = 0; xi < G_graphics_state->imgui_debug_texture_width; xi += 8)
                    {
                        __m256 x = _mm256_add_ps(x_base, x_accum);
                        
                        __m256 n = sample_terrain(x, y, z);
                        
                        f32 n_array[8];
                        _mm256_storeu_ps(n_array, n);
                        for(u32 pi = 0; pi < 8; pi++)
                        {
                            f32 pixel_noise = n_array[pi];
                            u32 c = pixel_noise < 0 ? 0xFF222222: 0xFF888888;
                            G_graphics_state->imgui_debug_texture_data[yi*G_graphics_state->imgui_debug_texture_width + xi + pi] = c;
                        }
                        
                        x_accum = _mm256_add_ps(x_accum, _mm256_set1_ps(8.0f));
                    }
                    y = _mm256_add_ps(y, _mm256_set1_ps(1.0f));
                }
                
                glBindTexture(GL_TEXTURE_2D, G_graphics_state->imgui_debug_texture);
                glTexSubImage2D(
                                GL_TEXTURE_2D,
                                0,
                                0,
                                0,
                                G_graphics_state->imgui_debug_texture_width,
                                G_graphics_state->imgui_debug_texture_height,
                                GL_RGBA,
                                GL_UNSIGNED_BYTE,
                                G_graphics_state->imgui_debug_texture_data);
                
                ImGui::Image((void*)(intptr_t)G_graphics_state->imgui_debug_texture,
                             ImVec2(f32(G_graphics_state->imgui_debug_texture_width), f32(G_graphics_state->imgui_debug_texture_width)));
                glBindTexture(GL_TEXTURE_2D, 0);
                
                ImGui::EndTabItem();
            }
            
            if(imgui_tab_bar && ImGui::BeginTabItem("Profile"))
            {
                for(u32 i = 0; i < G_profile_data.num; i++)
                {
                    ImGui::Text("%s: %i us", G_profile_data.name[i], G_profile_data.time[i]);
                }
                ImGui::EndTabItem();
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
                ImGui::Text(G_graphics_state->debug_info_log);
                
                G_profile_data.num = 0;
            }
            
            if(imgui_tab_bar)
            {
                ImGui::EndTabBar();
            }
            
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            
            glfwSwapBuffers(G_graphics_state->window);
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


