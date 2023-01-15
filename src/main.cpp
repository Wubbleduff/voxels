
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


// NOTE: This value must match the shader.
static const u32 SHADER_BUFFER_WIDTH = 1024*1024;
static const u32 NUM_VOXEL_DATA_FIELDS = 4;

static const u32 MAX_VOXELS = 10*1024*1024;
struct VoxelRenderData
{
    u32 num;
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

struct VoxelData
{
    u32 num;
    // | xxxx yyyy zzz |
    s32 pos[MAX_VOXELS*3];
    u32 color[MAX_VOXELS];
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
    const Camera* __restrict cam
)
{
    v3 plane_normals[6];
    v3 plane_points[6];
    get_frustum_planes(plane_normals, plane_points, cam);

    __m256 p_dot_n[] =
    {
        _mm256_set1_ps(dot(plane_points[0], plane_normals[0])),
        _mm256_set1_ps(dot(plane_points[1], plane_normals[1])),
        _mm256_set1_ps(dot(plane_points[2], plane_normals[2])),
        _mm256_set1_ps(dot(plane_points[3], plane_normals[3])),
        _mm256_set1_ps(dot(plane_points[4], plane_normals[4])),
        _mm256_set1_ps(dot(plane_points[5], plane_normals[5]))
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
            __m256 nx = _mm256_broadcast_ss(&plane_normals[ni].x);
            __m256 ny = _mm256_broadcast_ss(&plane_normals[ni].y);
            __m256 nz = _mm256_broadcast_ss(&plane_normals[ni].z);
            __m256 d = _mm256_fmadd_ps(vx, nx, _mm256_fmadd_ps(vy, ny, _mm256_mul_ps(vz, nz)));

            __m256 dist = _mm256_sub_ps(d, p_dot_n[ni]);
            __m256 cmp = _mm256_cmp_ps(dist, _mm256_set1_ps(1.0f), _CMP_LT_OQ);
            mask = _mm256_and_ps(mask, cmp);
        }

        u32 pack_index = _mm256_movemask_ps(mask);
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*0 + num_voxels,   left_pack(vx, pack_index));
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*1 + num_voxels,   left_pack(vy, pack_index));
        _mm256_storeu_ps(out_voxel_pos + MAX_VOXELS*2 + num_voxels,   left_pack(vz, pack_index));
        _mm256_storeu_si256((__m256i*)(out_voxel_color + num_voxels), left_pack(vc, pack_index));
        num_voxels += _mm_popcnt_u32(pack_index);
    }
    *out_num_voxels = num_voxels;
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

// Right now I'm just defining "chunk" as a MxNxO group of voxels
static void terrain_chunk(VoxelData *out, u32 chunk_x, u32 chunk_z, u32 width, u32 depth, f32 noise_scale_x, f32 noise_scale_y, f32 noise_scale_z)
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
    for(s32 z = 0; z < depth; z++)
    {
        for(s32 x = 0; x < width; x++)
        {
            s32 X = x + chunk_x;
            s32 Z = z + chunk_z;

            s32 Y = (s32)noise(X, Z);

            // Check for holes in the ground
            {
                s32 y0 = (s32)noise(X - 1, Z);
                s32 y1 = (s32)noise(X + 1, Z);
                s32 y2 = (s32)noise(X, Z - 1);
                s32 y3 = (s32)noise(X, Z + 1);
                s32 max_d = max(Y - y0, max(Y - y1, max(Y - y2, Y - y3)));
                for(s32 i = 0; i < max_d - 1; i++)
                {
                    s32 yi = Y - 1 - i;
                    out->pos[MAX_VOXELS*0 + out->num] = X;
                    out->pos[MAX_VOXELS*1 + out->num] = yi;
                    out->pos[MAX_VOXELS*2 + out->num] = Z;
                    out->color[out->num] = pick_color(yi);
                    out->num++;
                    if(out->num >= MAX_VOXELS) return;
                }
            }

            out->pos[MAX_VOXELS*0 + out->num] = X;
            out->pos[MAX_VOXELS*1 + out->num] = Y;
            out->pos[MAX_VOXELS*2 + out->num] = Z;
            out->color[out->num] = pick_color(Y);
            out->num++;
            if(out->num >= MAX_VOXELS) return;
        }
    }
}

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

    VoxelData * __restrict all_voxel_data = new VoxelData;
    VoxelRenderData * __restrict voxel_render_data = new VoxelRenderData;

    static f32 noise_scale_x = 0.008f;
    static f32 noise_scale_y = 20.0f;
    static f32 noise_scale_z = 0.008f;
    static const u32 WIDTH = 1024;
    static const u32 DEPTH = 1024;
    all_voxel_data->num = 0;
    for(s32 z = -1; z < 2; z++)
    {
        for(s32 x = -1; x < 2; x++)
        {
            terrain_chunk(all_voxel_data, x*WIDTH, z*DEPTH, WIDTH, DEPTH, noise_scale_x, noise_scale_y, noise_scale_z);
            if(all_voxel_data->num >= MAX_VOXELS) break;
        }
        if(all_voxel_data->num >= MAX_VOXELS) break;
    }

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
                static f32 camera_speed = 10.0f;
                ImGui::DragFloat("camera speed", &camera_speed);
                current_cam->pos += normalize(camera_vel) * TIME_STEP * camera_speed;
            }

            s32 camera_x = (s32)graphics_state->cam.pos.x;
            s32 camera_z = (s32)graphics_state->cam.pos.z;
            s32 camera_chunk_x = camera_x & ~0b0001'1111'1111;
            s32 camera_chunk_z = camera_z & ~0b0001'1111'1111;
            bool should_gen = false;
            static s32 last_camera_chunk_x = 0xFFFFFFFF;
            static s32 last_camera_chunk_z = 0xFFFFFFFF;
            /*
            if((camera_chunk_x != last_camera_chunk_x) || (camera_chunk_z != last_camera_chunk_z))
            {
                should_gen = true;
            }
            */
            if(should_gen)
            {
                all_voxel_data->num = 0;
                for(s32 z = -1; z <= 1; z++)
                {
                    for(s32 x = -1; x <= 1; x++)
                    {
                        terrain_chunk(all_voxel_data,
                                camera_chunk_x + x*WIDTH,
                                camera_chunk_z + z*DEPTH,
                                WIDTH, DEPTH, noise_scale_x, noise_scale_y, noise_scale_z);
                        if(all_voxel_data->num >= MAX_VOXELS) break;
                    }
                    if(all_voxel_data->num >= MAX_VOXELS) break;
                }
            }
            last_camera_chunk_x = camera_chunk_x;
            last_camera_chunk_z = camera_chunk_z;

            // Prep render data
            {
                TIME_SCOPE("Prep render data");

                ImGui::Text("prep num_voxels %i", all_voxel_data->num);

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
                const mat4 view_mat = view_m_world(&graphics_state->cam);
                float camera_xform[] =
                {
                    view_mat[0][0], view_mat[0][1], view_mat[0][2],
                    view_mat[1][0], view_mat[1][1], view_mat[1][2],
                    view_mat[2][0], view_mat[2][1], view_mat[2][2],
                    graphics_state->cam.pos.x, graphics_state->cam.pos.y, graphics_state->cam.pos.z
                };
                camera_cull(
                    &voxel_render_data->num,
                    voxel_render_data->pos,
                    voxel_render_data->color,
                    all_voxel_data->num,
                    all_voxel_data->pos,
                    all_voxel_data->color,
                    &graphics_state->cam
                );
#endif
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

