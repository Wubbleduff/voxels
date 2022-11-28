
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

struct VoxelRenderData
{
    static const u32 MAX_VOXELS = 10*1024*1024;
    u32 num;
    f32 x[MAX_VOXELS];
    f32 y[MAX_VOXELS];
    f32 z[MAX_VOXELS];
    u32 color[MAX_VOXELS];

    static const u32 BATCH_SIZE = 10*1024;
};

struct GraphicsState
{
    GLFWwindow *window;
    char debug_info_log[1024];
    
    f32 screen_aspect_ratio;

    GLuint batch_voxel_shader_program;
    GLuint batch_voxel_vao;
    GLuint batch_voxel_vbo;
    GLuint batch_voxel_ebo;
    GLuint batch_voxel_ssbo;

    v3 camera_pos;
    f32 camera_x_rot;
    f32 camera_y_rot;
    //f32 camera_width;
    f32 camera_fov;
    f32 view_dist;
};
static GraphicsState *graphics_state = nullptr;

struct InputState
{
    //bool keys[256];
    s32 mouse_screen_x;
    s32 mouse_screen_y;
};
static InputState *input_state = nullptr;

struct VoxelData
{
    static const u32 MAX_VOXELS = 10*1024*1024;
    u32 num;
    s32 x[MAX_VOXELS];
    s32 y[MAX_VOXELS];
    s32 z[MAX_VOXELS];

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

GLuint make_shader(const char *shader_path)
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
    
    u32 vert_shader;
    vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vert_source, NULL);
    glCompileShader(vert_shader);
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(vert_shader, info_log_size, NULL, graphics_state->debug_info_log);
        //fprintf(stderr, "VERTEX SHADER ERROR : %s\n%s\n", shader_path, info_log);
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
        //fprintf(stderr, "VERTEX SHADER ERROR : %s\n%s\n", shader_path, info_log);
        return -1;
    }
    
    
    u32 geom_shader;
    if(geom_source != nullptr)
    {
        geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        glShaderSource(geom_shader, 1, &geom_source, NULL);
        glCompileShader(geom_shader);
        glGetShaderiv(geom_shader, GL_COMPILE_STATUS, &success);
        if(!success)
        {
            glGetShaderInfoLog(geom_shader, info_log_size, NULL, graphics_state->debug_info_log);
            //fprintf(stderr, "VERTEX SHADER ERROR : %s\n%s\n", shader_path, info_log);
            return -1;
        }
    }
    
    check_gl_errors("compiling shaders");
    
    GLuint program = glCreateProgram();
    check_gl_errors("making program");
    
    glAttachShader(program, vert_shader);
    glAttachShader(program, frag_shader);
    if(geom_source != nullptr) glAttachShader(program, geom_shader);
    glLinkProgram(program);
    check_gl_errors("linking program");
    
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success)
    {
        glGetProgramInfoLog(program, info_log_size, NULL, graphics_state->debug_info_log);
        //fprintf(stderr, "SHADER LINKING ERROR : %s\n%s\n", shader_path, info_log);
        return -1;
    }
    
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader); 
    check_gl_errors("deleting shaders");
    
    return program;
}

static void reshape(GLFWwindow *window, s32 width, s32 height)
{
    glViewport(0, 0, (GLint)width, (GLint)height);
    graphics_state->screen_aspect_ratio = (f32)width / height;
}

// TODO speed
mat4 view_m_world()
{

    mat4 y_rot = make_y_axis_rotation_matrix(graphics_state->camera_y_rot);
    mat4 x_rot = make_x_axis_rotation_matrix(graphics_state->camera_x_rot);

    //return make_translation_matrix(-graphics_state->camera_pos) * x_rot * y_rot;
    return x_rot * y_rot * make_translation_matrix(-graphics_state->camera_pos);
}
mat4 clip_m_view()
{
    static f32 n = 1.0f;
    static f32 f = graphics_state->view_dist;
    f32 fov = deg_to_rad(graphics_state->camera_fov);
    f32 r = -(f + n) / (f - n);
    f32 s = -(2.0f * n * f) / (f - n);
    mat4 result =
    {
         (f32)(1.0f / tanf(fov / 2.0f)) / graphics_state->screen_aspect_ratio, 0.0f, 0.0f, 0.0f,
         0.0f, 1.0f / tanf(fov / 2.0f), 0.0f, 0.0f,
         0.0f, 0.0f, r, s,
         0.0f, 0.0f, -1.0f, 0.0f
    };
    
    return result;
}

static void terrain(VoxelData *out, u32 width, u32 depth, float noise_scale_x, float noise_scale_y, float noise_scale_z)
{
    out->num = 0;

    auto noise = [noise_scale_x, noise_scale_y, noise_scale_z](s32 x_, s32 z_)
    {
        //return (perlin_noise(x_*noise_scale_x, z_*noise_scale_z)*0.5f + 0.5f)*noise_scale_y;
        float n = (perlin_noise(x_*noise_scale_x, z_*noise_scale_z)*0.5f + 0.5f)*noise_scale_y;
        n += (perlin_noise(x_*noise_scale_x*0.1f, z_*noise_scale_z*0.1f)*0.5f + 0.5f)*noise_scale_y*10.0f;

        const static s32 s_x[] = {128, 512, -8*128,  7*128, 7*128, -8*128};
        const static s32 s_z[] = {128, 512, -8*128, -8*128, 7*128,  7*128};
        const static u32 num_s = sizeof(s_x) / sizeof(s_x[0]);
        s32 min_d = 1000000000; // TODO
        for(u32 i = 0; i < num_s; i++)
        {
            s32 d_x = x_ - s_x[i];
            s32 d_z = z_ - s_z[i];
            s32 d = (s32)sqrtf(d_x*d_x + d_z*d_z);
            if(d < min_d)
            {
                min_d = d;
            }
        }
        n += min_d*0.3f;
        
        return n;
    };
    auto pick_color = [](s32 y)
    {
        if(y < 300) return (u32)0x424242FF;
        else if(y < 480) return (u32)0x519916FF;
        else return (u32)0xCCCCCCFF;
    };
    for(s32 chunk_z = -8; chunk_z < 8; chunk_z++)
    {
        for(s32 chunk_x = -8; chunk_x < 8; chunk_x++)
        {
            for(s32 z = 0; z < depth; z++)
            {
                for(s32 x = 0; x < width; x++)
                {
                    s32 X = x + (chunk_x)*width;
                    s32 Z = z + (chunk_z)*depth;

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
                            out->x[out->num] = X;
                            out->y[out->num] = yi;
                            out->z[out->num] = Z;
                            out->color[out->num] = pick_color(yi);
                            out->num++;
                            if(out->num >= VoxelData::MAX_VOXELS) return;
                        }
                    }

                    out->x[out->num] = X;
                    out->y[out->num] = Y;
                    out->z[out->num] = Z;
                    out->color[out->num] = pick_color(Y);
                    out->num++;
                    if(out->num >= VoxelData::MAX_VOXELS) return;
                }
            }
        }
    }
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

        graphics_state->camera_pos = v3(0.0f, 32.0f, 0.0f);
        graphics_state->camera_x_rot = 0.5f;
        graphics_state->camera_fov = 90.0f;
        graphics_state->view_dist = 5000.0f;

        glGenVertexArrays(1, &graphics_state->batch_voxel_vao);
        glBindVertexArray(graphics_state->batch_voxel_vao);
        check_gl_errors("vao");

        // https://www.khronos.org/opengl/wiki/Shader_Storage_Buffer_Object
        glGenBuffers(1, &graphics_state->batch_voxel_vbo);
        check_gl_errors("vbo");
        static f32 voxel_vertices[] =
        {
             0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,
             0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f, 
            -0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f,
        };
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
        glBufferData(GL_ARRAY_BUFFER, sizeof(voxel_vertices), voxel_vertices, GL_STATIC_DRAW);
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

        graphics_state->batch_voxel_shader_program = make_shader("assets/shaders/batch_voxel.shader");
        {
            FILE *file = fopen("log.txt", "wt");
            fprintf(file, graphics_state->debug_info_log);
            fclose(file);
        }
        //assert(graphics_state->batch_voxel_shader_program != -1);

        const u32 voxel_buffer_size = SHADER_BUFFER_WIDTH * NUM_VOXEL_DATA_FIELDS * 4;
        glGenBuffers(1, &graphics_state->batch_voxel_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, graphics_state->batch_voxel_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, voxel_buffer_size, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, graphics_state->batch_voxel_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        check_gl_errors("ssbo");
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
    static f32 noise_scale_y = 48.0f;
    static f32 noise_scale_z = 0.008f;
    static const u32 WIDTH = 128;
    static const u32 DEPTH = 128;
    terrain(all_voxel_data, WIDTH, DEPTH, noise_scale_x, noise_scale_y, noise_scale_z);

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
                // glClearColor(0.01f, 0.0f, 0.1f, 0.0f); // Deep blue
                glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                ImGui_ImplOpenGL3_NewFrame();
                ImGui_ImplGlfw_NewFrame();
                ImGui::NewFrame();
            }

            running = !glfwGetKey(graphics_state->window, GLFW_KEY_ESCAPE) && !glfwWindowShouldClose(graphics_state->window);

            {
                static f32 last_mouse_x, last_mouse_y;
                double xpos, ypos;
                glfwGetCursorPos(graphics_state->window, &xpos, &ypos);
                if(glfwGetKey(graphics_state->window, GLFW_KEY_SPACE))
                {
                    v2 mouse_delta = v2((f32)xpos - last_mouse_x, (f32)ypos - last_mouse_y);
                    graphics_state->camera_y_rot += mouse_delta.x * TIME_STEP * 0.2f;
                    graphics_state->camera_x_rot += mouse_delta.y * TIME_STEP * 0.2f;
                }
                last_mouse_x = (f32)xpos;
                last_mouse_y = (f32)ypos;

                mat4 y_rot = make_y_axis_rotation_matrix(graphics_state->camera_y_rot);
                mat4 x_rot = make_x_axis_rotation_matrix(graphics_state->camera_x_rot);
                mat4 camera_t = x_rot * y_rot;
                v3 camera_right   = v3(camera_t[0][0], camera_t[0][1], camera_t[0][2]);
                v3 camera_forward = v3(camera_t[2][0], camera_t[2][1], camera_t[2][2]);

                v3 camera_vel = v3();
                camera_vel += -camera_right   * glfwGetKey(graphics_state->window, GLFW_KEY_A);
                camera_vel +=  camera_right   * glfwGetKey(graphics_state->window, GLFW_KEY_D);
                camera_vel +=  camera_forward * glfwGetKey(graphics_state->window, GLFW_KEY_S);
                camera_vel += -camera_forward * glfwGetKey(graphics_state->window, GLFW_KEY_W);
                static f32 camera_speed = 100.0f;
                ImGui::DragFloat("camera speed", &camera_speed);
                graphics_state->camera_pos += normalize(camera_vel) * TIME_STEP * camera_speed;

                //graphics_state->camera_pos.y = (s32)(noise(graphics_state->camera_pos.x, graphics_state->camera_pos.z) + 2.0f);
            }


            // TODO cleanup
            //static v3 ground_color = v3(93.0f, 85.0f, 60.0f);
            static float view_dist = 5000.0f;
            {
                ImGui::DragFloat("noise_scale_x", &noise_scale_x);
                ImGui::DragFloat("noise_scale_y", &noise_scale_y);
                ImGui::DragFloat("noise_scale_z", &noise_scale_z);
                if(ImGui::Button("gen"))
                {
                    terrain(all_voxel_data, WIDTH, DEPTH, noise_scale_x, noise_scale_y, noise_scale_z);
                }
                ImGui::DragFloat("view dist", &view_dist);
            }

            // Prep render data
            {
                TIME_SCOPE("Prep render data");
                mat4 y_rot = make_y_axis_rotation_matrix(graphics_state->camera_y_rot);
                mat4 x_rot = make_x_axis_rotation_matrix(graphics_state->camera_x_rot);
                mat4 camera_t = x_rot * y_rot;
                v3 camera_forward = v3(camera_t[2][0], camera_t[2][1], camera_t[2][2]);

                const float view_dist_sq = view_dist*view_dist;
                __m256 cf_x = _mm256_set1_ps(camera_forward.x);
                __m256 cf_y = _mm256_set1_ps(camera_forward.y);
                __m256 cf_z = _mm256_set1_ps(camera_forward.z);
                voxel_render_data->num = 0;
                for(u32 i = 0; i < all_voxel_data->num; i += 8)
                {
                    __m256 v_x = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(all_voxel_data->x + i)));
                    __m256 v_y = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(all_voxel_data->y + i)));
                    __m256 v_z = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i*)(all_voxel_data->z + i)));
                    __m256i v_color = _mm256_loadu_si256((__m256i*)(all_voxel_data->color + i));
                    __m256 c_x = _mm256_set1_ps(graphics_state->camera_pos.x);
                    __m256 c_y = _mm256_set1_ps(graphics_state->camera_pos.y);
                    __m256 c_z = _mm256_set1_ps(graphics_state->camera_pos.z);
                    __m256 mask = _mm256_cvtepi32_ps(_mm256_set1_epi32(0xFFFFFFFF));

                    __m256 d_x = _mm256_sub_ps(v_x, c_x);
                    __m256 d_y = _mm256_sub_ps(v_y, c_y);
                    __m256 d_z = _mm256_sub_ps(v_z, c_z);

                    __m256 len_sq = _mm256_fmadd_ps(d_x,d_x, _mm256_fmadd_ps(d_y,d_y, _mm256_mul_ps(d_z,d_z)));
                    mask = _mm256_and_ps(mask, _mm256_cmp_ps(len_sq, _mm256_set1_ps(view_dist_sq), _CMP_NGT_UQ));

                    __m256 camera_align = _mm256_fmadd_ps(d_x,cf_x, _mm256_fmadd_ps(d_y,cf_y, _mm256_mul_ps(d_z,cf_z)));
                    mask = _mm256_and_ps(mask, _mm256_cmp_ps(camera_align, _mm256_set1_ps(0.0f), _CMP_NGT_UQ));

                    s32 mask32 = _mm256_movemask_ps(mask);
                    u64 expanded_mask = _pdep_u64(mask32, 0x0101010101010101);
                    expanded_mask *= 0xFF;
                    u64 wanted_indices = _pext_u64(0x0706050403020100, expanded_mask);
                    __m256i shufmask = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(wanted_indices));

                    __m256 packed_x = _mm256_permutevar8x32_ps(v_x, shufmask);
                    __m256 packed_y = _mm256_permutevar8x32_ps(v_y, shufmask);
                    __m256 packed_z = _mm256_permutevar8x32_ps(v_z, shufmask);
                    __m256i packed_color = _mm256_permutevar8x32_epi32(v_color, shufmask);
                    _mm256_storeu_ps(&(voxel_render_data->x[voxel_render_data->num]), packed_x);
                    _mm256_storeu_ps(&(voxel_render_data->y[voxel_render_data->num]), packed_y);
                    _mm256_storeu_ps(&(voxel_render_data->z[voxel_render_data->num]), packed_z);
                    _mm256_storeu_si256((__m256i*)&(voxel_render_data->color[voxel_render_data->num]), packed_color);
                    voxel_render_data->num += _mm_popcnt_u64(mask32);
                }
            }

            // Draw
            {
                TIME_SCOPE("Draw");

                u32 num_batches_drawn = 0;
                for(u32 batch_i = 0; batch_i < voxel_render_data->num; batch_i += VoxelRenderData::BATCH_SIZE)
                {
                    float *batch_x = voxel_render_data->x + batch_i;
                    float *batch_y = voxel_render_data->y + batch_i;
                    float *batch_z = voxel_render_data->z + batch_i;
                    u32 *batch_color = voxel_render_data->color + batch_i;
                    u32 batch_size = min(voxel_render_data->num - batch_i, VoxelRenderData::BATCH_SIZE);

                    assert(VoxelRenderData::BATCH_SIZE <= SHADER_BUFFER_WIDTH);
                    glUseProgram(graphics_state->batch_voxel_shader_program);
                    check_gl_errors("use program");

                    {
                        mat4 m = view_m_world();
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "m_view");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }
                    {
                        mat4 m = clip_m_view();
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "m_proj");
                        glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
                        if(loc == -1) assert(false);
                    }
                    /*
                    {
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "u_view_dist");
                        glUniform1f(loc, view_dist);
                        if(loc == -1) assert(false);
                    }
                    */
                    /*
                    {
                        GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "u_color");
                        glUniform3f(loc, ground_color.x / 360.0f, ground_color.y / 100.0f, ground_color.z / 100.0f);
                        if(loc == -1) assert(false);
                    }
                    */

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

