
#include <windows.h>
#include <stdio.h> // fopen

#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "common.h"
#include "game_math.h"

// TODO  move
static const u32 chunk_w = 128;
static const u32 chunk_h = 32;
static const u32 chunk_d = 128;
static const u32 MAX_VOXELS = chunk_w*chunk_h*chunk_d;
static const u32 NUM_VOXEL_DATA_FIELDS = 3;

static const u32 SHADER_BUFFER_WIDTH = 65536;

struct GraphicsState
{
    GLFWwindow *window;
    
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
};
static GraphicsState *graphics_state = nullptr;

struct InputState
{
    //bool keys[256];
    s32 mouse_screen_x;
    s32 mouse_screen_y;
};
static InputState *input_state = nullptr;

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
    s32  success;
    char info_log[512];
    
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
        glGetShaderInfoLog(vert_shader, 512, NULL, info_log);
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
        glGetShaderInfoLog(frag_shader, 512, NULL, info_log);
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
            glGetShaderInfoLog(geom_shader, 512, NULL, info_log);
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
        glGetProgramInfoLog(program, 512, NULL, info_log);
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
mat4 ndc_m_world()
{
#if 0
    mat4 ndc_m_view =
    {
        2.0f / graphics_state->camera_width, 0.0f, 0.0f, 0.0f,
        0.0f, (2.0f / graphics_state->camera_width) * graphics_state->screen_aspect_ratio, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
#else

    static f32 n = 1.0f;
    static f32 f = 100.0f;
    f32 fov = deg_to_rad(graphics_state->camera_fov);
    f32 r = -(f + n) / (f - n);
    f32 s = -(2.0f * n * f) / (f - n);
    mat4 ndc_m_view =
    {
         (f32)(1.0f / tanf(fov / 2.0f)) / graphics_state->screen_aspect_ratio, 0.0f, 0.0f, 0.0f,
         0.0f, 1.0f / tanf(fov / 2.0f), 0.0f, 0.0f,
         0.0f, 0.0f, r, s,
         0.0f, 0.0f, -1.0f, 0.0f
    };
#endif
    
    
    return ndc_m_view * view_m_world();
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

        graphics_state->camera_pos = v3(0.0f, 2.0f, 0.0f);
        graphics_state->camera_x_rot = 0.5f;
        graphics_state->camera_fov = 90.0f;

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
        assert(graphics_state->batch_voxel_shader_program != -1);

        const u32 voxel_buffer_size = SHADER_BUFFER_WIDTH * NUM_VOXEL_DATA_FIELDS * sizeof(f32);
        //assert(MAX_VOXELS * NUM_VOXEL_DATA_FIELDS * sizeof(f32) < voxel_buffer_size);
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
                glClearColor(0.01f, 0.0f, 0.1f, 0.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                ImGui_ImplOpenGL3_NewFrame();
                ImGui_ImplGlfw_NewFrame();
                ImGui::NewFrame();
            }

            running = !glfwGetKey(graphics_state->window, GLFW_KEY_ESCAPE) && !glfwWindowShouldClose(graphics_state->window);

            {
                static float last_mouse_x, last_mouse_y;
                double xpos, ypos;
                glfwGetCursorPos(graphics_state->window, &xpos, &ypos);
                if(glfwGetKey(graphics_state->window, GLFW_KEY_SPACE))
                {
                    v2 mouse_delta = v2((float)xpos - last_mouse_x, (float)ypos - last_mouse_y);
                    graphics_state->camera_y_rot += mouse_delta.x * TIME_STEP * 0.2f;
                    graphics_state->camera_x_rot += mouse_delta.y * TIME_STEP * 0.2f;
                }
                last_mouse_x = (float)xpos;
                last_mouse_y = (float)ypos;

                mat4 y_rot = make_y_axis_rotation_matrix(graphics_state->camera_y_rot);
                mat4 x_rot = make_x_axis_rotation_matrix(graphics_state->camera_x_rot);
                mat4 asdf = x_rot * y_rot;

                //v3 camera_right   = v3(1.0f, 0.0f, 0.0f);
                //v3 camera_forward = v3(0.0f, 0.0f, 1.0f);
                v3 camera_right   = v3(asdf[0][0], asdf[0][1], asdf[0][2]);
                v3 camera_forward = v3(asdf[2][0], asdf[2][1], asdf[2][2]);

                v3 camera_vel = v3();
                camera_vel += -camera_right   * glfwGetKey(graphics_state->window, GLFW_KEY_A);
                camera_vel +=  camera_right   * glfwGetKey(graphics_state->window, GLFW_KEY_D);
                camera_vel +=  camera_forward * glfwGetKey(graphics_state->window, GLFW_KEY_S);
                camera_vel += -camera_forward * glfwGetKey(graphics_state->window, GLFW_KEY_W);
                graphics_state->camera_pos += normalize(camera_vel) * TIME_STEP * 30.0f;
            }

            static bool gen = true;
            static f32 ox = 0.025f;
            static f32 oy = 0.025f;
            static f32 scale = 1.0f;
            gen = ImGui::DragFloat("x", &ox, 0.001f) || gen;
            gen = ImGui::DragFloat("y", &oy, 0.001f) || gen;
            gen = ImGui::DragFloat("scale", &scale) || gen;
            ldif(gen)
            {
                u32 num_voxels = 0;
                static f32 data[NUM_VOXEL_DATA_FIELDS][MAX_VOXELS] = {};

                glUseProgram(graphics_state->batch_voxel_shader_program);
                check_gl_errors("use program");

                mat4 vp = ndc_m_world();
                GLint loc = glGetUniformLocation(graphics_state->batch_voxel_shader_program, "vp");
                glUniformMatrix4fv(loc, 1, true, &(vp[0][0]));
                if(loc == -1) assert(false);

                glBindVertexArray(graphics_state->batch_voxel_vao);
                glBindBuffer(GL_ARRAY_BUFFER, graphics_state->batch_voxel_vbo);

                for(u32 chunk_i = 0; chunk_i < 2; chunk_i++)
                {
                    for(s32 y = 0; y < chunk_h; y++)
                    {
                        for(s32 z = 0; z < chunk_d; z++)
                        {
                            for(s32 x = 0; x < chunk_w; x++)
                            {
                                s32 X = x + chunk_i*chunk_w;
                                //s32 X = x;
                                f32 height = max(perlin_noise(X*ox, z*oy)*scale + 2.0f, 2.0f);

                                s32 x0 = X - 1;
                                s32 x1 = X + 1;
                                s32 z0 = z - 1;
                                s32 z1 = z + 1;

                                f32 height_x0 = max(perlin_noise(x0*ox, z*oy)*scale + 2.0f, 2.0f);
                                f32 height_x1 = max(perlin_noise(x1*ox, z*oy)*scale + 2.0f, 2.0f);
                                f32 height_z0 = max(perlin_noise(X*ox, z0*oy)*scale + 2.0f, 2.0f);
                                f32 height_z1 = max(perlin_noise(X*ox, z1*oy)*scale + 2.0f, 2.0f);

                                float Y = (f32)y;
                                bool surrounded =
                                    Y < height_x0 && Y < height_x1 &&
                                    Y + 1 < height && Y - 1 < height &&
                                    Y < height_z0 && Y < height_z1;

                                if((f32)y < height && !surrounded)
                                {
                                    data[0][num_voxels] = (f32)X;
                                    data[1][num_voxels] = (f32)y;
                                    data[2][num_voxels] = (f32)z;
                                    num_voxels++;
                                }
                            }
                        }
                    }

                    u32 offset_x = 0*SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_y = 1*SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 offset_z = 2*SHADER_BUFFER_WIDTH*sizeof(f32);
                    u32 size = sizeof(f32) * num_voxels;
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, graphics_state->batch_voxel_ssbo);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_x, size, data[0]);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_y, size, data[1]);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset_z, size, data[2]);
                    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, graphics_state->batch_voxel_ssbo);

                    ImGui::Text("num_voxels %i", num_voxels);

                    glDrawElementsInstanced(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0, num_voxels);
                    check_gl_errors("draw");
                }

            }
            gen = false;

            {
                f64 frame_ms = (glfwGetTime() - current_time) * 1000.0;
                static f64 frame_times[32] = {};
                static u32 num_frame_times = 0;
                frame_times[num_frame_times++] = frame_ms;
                if(num_frame_times >= 32) num_frame_times = 0;
                f64 frame_time_sum = 0.0;
                for(u32 i = 0; i < 32; i++) frame_time_sum += frame_times[i];
                ImGui::Text("frame time: %f ms", frame_time_sum / 32.0);
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

