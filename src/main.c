
#include "common.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <GL/GL.h>
#include "wglext.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Globals

HANDLE g_log_file;
struct OpenGLState* g_opengl_state;
struct InputState* g_input_state;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

typedef struct
{
    // Last byte should always be 0 so we can trivially cast to C-string.
    _Alignas(32) char buf[32];
} StringBuf32;

static inline const char* StringBuf32_to_cstr(const StringBuf32* s)
{
    return s->buf;
}

static inline StringBuf32 u32_to_StringBuf32(u32_m n)
{
    StringBuf32 result;
    memset(result.buf, 0, sizeof(result));
    u64_m size = 0;
    while(n && size < 31)
    {
        u32 digit = n % 10;
        result.buf[size++] = (char)('0' + digit);
        n /= 10;
    }
    for(u64_m i = 0; i < size / 2; i++)
    {
        char tmp = result.buf[i];
        result.buf[i] = result.buf[size - i - 1];
        result.buf[size - i - 1] = tmp;
    }
    return result;
}


static inline u64 cstr_len(const char* s)
{
    u64_m n = 0;
    while(*s++) n++;
    return n;
}

static inline void assert_fn(const char* file, int line, s32 c, const char* msg)
{
    if(!c)
    {
        const StringBuf32 line_str = u32_to_StringBuf32(line);
        const char* line_cstr = StringBuf32_to_cstr(&line_str);
        WriteFile(g_log_file, file, (int)cstr_len(file), NULL, NULL);
        WriteFile(g_log_file, ":", 1, NULL, NULL);
        WriteFile(g_log_file, line_cstr, (int)cstr_len(line_cstr), NULL, NULL);
        WriteFile(g_log_file, "  ", 1, NULL, NULL);
        WriteFile(g_log_file, msg, (int)cstr_len(msg), NULL, NULL);
        FlushFileBuffers(g_log_file);
        DebugBreak();
    }
}
#define ASSERT(c, msg) assert_fn(__FILE__, __LINE__, (c), (msg))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL / WGL

// Definitions copied from glad library.
typedef signed long long int GLsizeiptr;
typedef signed long long int GLintptr;
typedef char GLchar;
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_ARRAY_BUFFER 0x8892
#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82

// https://www.khronos.org/opengl/wiki/Load_OpenGL_Functions
void *load_gl_fn(HMODULE opengl32_dll_module, const char *name)
{
    void* p = (void*)wglGetProcAddress(name);
    if(p == 0 || (p == (void*)0x1) || (p == (void*)0x2) || (p == (void*)0x3) || (p == (void*)-1) )
    {
        p = (void *)GetProcAddress(opengl32_dll_module, name);
    }

    ASSERT(p != NULL, "Could not load OpenGL function.");
    ASSERT(p != (void*)1, "Could not load OpenGL function.");
    ASSERT(p != (void*)2, "Could not load OpenGL function.");
    ASSERT(p != (void*)-1, "Could not load OpenGL function.");

    return p;
}

typedef HGLRC (*fnptr_wglCreateContextAttribsARB)(HDC hDC, HGLRC hshareContext, const int *attribList);
typedef void (*fnptr_glGenVertexArrays)(GLsizei n, GLuint *arrays);
typedef GLenum (*fnptr_glGetError)();
typedef void (*fnptr_glBindVertexArray)(GLuint array);
typedef void (*fnptr_glGenBuffers)(GLsizei n, GLuint * buffers);
typedef void (*fnptr_glBindBuffer)(GLenum target, GLuint buffer);
typedef void (*fnptr_glBufferData)(GLenum target, GLsizeiptr size, const void * data, GLenum usage);
typedef void (*fnptr_glVertexAttribPointer)(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void * pointer);
typedef void (*fnptr_glEnableVertexAttribArray)(GLuint index);
typedef GLuint (*fnptr_glCreateShader)(GLenum shaderType);
typedef void (*fnptr_glShaderSource)(GLuint shader, GLsizei count, const GLchar **string, const GLint *length);
typedef void (*fnptr_glCompileShader)(GLuint shader);
typedef void (*fnptr_glGetShaderiv)(GLuint shader, GLenum pname, GLint *params);
typedef void (*fnptr_glGetShaderInfoLog)(GLuint shader, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
typedef GLuint (*fnptr_glCreateProgram)(void);
typedef void (*fnptr_glAttachShader)(GLuint program, GLuint shader);
typedef void (*fnptr_glLinkProgram)(GLuint program);
typedef void (*fnptr_glGetProgramiv)(GLuint program, GLenum pname, GLint *params);
typedef void (*fnptr_glGetProgramInfoLog)(GLuint program, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
typedef void (*fnptr_glDeleteShader)(GLuint shader);
typedef void (*fnptr_glUseProgram)(GLuint program);
typedef void (*fnptr_glBufferSubData)(GLenum target, GLintptr offset, GLsizeiptr size, const void * data);
typedef void (*fnptr_glDrawElementsInstanced)(GLenum mode, GLsizei count, GLenum type, const void * indices, GLsizei instancecount);

#define VERTEX_ARRAY_BYTES MB(1)
#define INDEX_ARRAY_BYTES MB(1)

struct OpenGLState
{
    HWND hwnd;
    HGLRC gl_context;
    GLuint last_gl_error;

    GLuint vertex_array_object;
    GLuint index_buffer_object;
    GLuint vertex_buffer_object_x;
    GLuint vertex_buffer_object_y;
    GLuint vertex_buffer_object_z;
    GLuint shader_program;

    fnptr_wglCreateContextAttribsARB wglCreateContextAttribsARB;
    fnptr_glGetError glGetError;
    fnptr_glGenVertexArrays glGenVertexArrays;
    fnptr_glBindVertexArray glBindVertexArray;
    fnptr_glGenBuffers glGenBuffers;
    fnptr_glBindBuffer glBindBuffer;
    fnptr_glBufferData glBufferData;
    fnptr_glVertexAttribPointer glVertexAttribPointer;
    fnptr_glEnableVertexAttribArray glEnableVertexAttribArray;
    fnptr_glCreateShader glCreateShader;
    fnptr_glShaderSource glShaderSource;
    fnptr_glCompileShader glCompileShader;
    fnptr_glGetShaderiv glGetShaderiv;
    fnptr_glGetShaderInfoLog glGetShaderInfoLog;
    fnptr_glCreateProgram glCreateProgram;
    fnptr_glAttachShader glAttachShader;
    fnptr_glLinkProgram glLinkProgram;
    fnptr_glGetProgramiv glGetProgramiv;
    fnptr_glGetProgramInfoLog glGetProgramInfoLog;
    fnptr_glDeleteShader glDeleteShader;
    fnptr_glUseProgram glUseProgram;
    fnptr_glBufferSubData glBufferSubData;
    fnptr_glDrawElementsInstanced glDrawElementsInstanced;
};

#define CALL_GL(fn, ...) \
do { \
    g_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = g_opengl_state->glGetError(); \
    g_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
} while(0)

#define CALL_GL_RET(out, type, fn, ...) \
do { \
    type r = g_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = g_opengl_state->glGetError(); \
    g_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
    *(out) = r; \
} while(0)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Input

#define MAX_KEYBOARD_KEYS 256
#define MAX_MOUSE_KEYS 32
struct InputState
{
    u32_m key[MAX_KEYBOARD_KEYS];
    u32_m mouse_key[MAX_MOUSE_KEYS];
};

enum KeyboardKey
{
    KB_ESCAPE = VK_ESCAPE,
};

static inline u32 is_key_down(enum KeyboardKey k)
{
    ASSERT(k >= 0, "is_key_down key code underflow");
    ASSERT(k < MAX_KEYBOARD_KEYS, "is_key_down key code overflow");
    return g_input_state->key[(u32)k];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



LRESULT CALLBACK WindowProc(HWND window, UINT message, WPARAM wParam, LPARAM lParam)
{
    LRESULT result = 0; 

    switch (message)
    {
        case WM_KEYDOWN: 
        {
            g_input_state->key[wParam] = 1;
            break;
        }

        case WM_KEYUP:
        {
            g_input_state->key[wParam] = 0;
            break;
        }

        case WM_LBUTTONDOWN:
        {
            g_input_state->mouse_key[0] = 1;
            break;
        }

        case WM_LBUTTONUP:
        {
            g_input_state->mouse_key[0] = 0;
            break;
        }

        case WM_RBUTTONDOWN:
        {
            g_input_state->mouse_key[1] = 1;
            break;
        }
        case WM_RBUTTONUP:
        {
            g_input_state->mouse_key[1] = 0;
            break;
        }

        case WM_DESTROY:
        {
            // https://learn.microsoft.com/en-us/windows/win32/learnwin32/closing-the-window?redirectedfrom=MSDN
            PostQuitMessage(0);
            return 0;
        }  

        default:
        {
            result = DefWindowProc(window, message, wParam, lParam);
            break;
        }
    }

    return result;
}

int WinMainCRTStartup()
{
    g_log_file = CreateFileA("log.txt", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    ASSERT(g_log_file != INVALID_HANDLE_VALUE, "Couldn't open log file for writing.");

    struct InputState input_state_storage;
    g_input_state = &input_state_storage;

    struct OpenGLState opengl_state_storage;
    g_opengl_state = &opengl_state_storage;

    {
        const HMODULE hinstance = GetModuleHandle(0);

        WNDCLASS window_class = {0};
        window_class.style = CS_HREDRAW|CS_VREDRAW|CS_OWNDC;
        window_class.lpfnWndProc = WindowProc;
        window_class.hInstance = hinstance;
        window_class.hCursor = LoadCursor(NULL, IDC_ARROW);
        window_class.lpszClassName = "Windows Program Class";
        ATOM register_class_result = RegisterClass(&window_class);
        ASSERT(register_class_result, "Register class failed");

        u64 window_width = GetSystemMetrics(SM_CXSCREEN) / 2;
        u64 window_height = GetSystemMetrics(SM_CYSCREEN) / 2;
        g_opengl_state->hwnd = CreateWindowEx(0,                     // Extended style
                window_class.lpszClassName,               // Class name
                "",                                       // Window name
                WS_OVERLAPPEDWINDOW | WS_VISIBLE,         // Style of the window
                0,                                        // Initial X position
                0,                                        // Initial Y position
                (int)window_width,                        // Initial width
                (int)window_height,                       // Initial height
                0,                                        // Handle to the window parent
                0,                                        // Handle to a menu
                hinstance,                                // Handle to an instance
                0);
        ASSERT(g_opengl_state->hwnd != NULL, "Failed to create a window");


        const HDC dc = GetDC(g_opengl_state->hwnd);
        ASSERT(dc != NULL, "GetDC failed.");

        HMODULE opengl32_dll_module = LoadLibraryA("opengl32.dll");
        ASSERT(opengl32_dll_module != NULL, "Could not load opengl32.dll");

        PIXELFORMATDESCRIPTOR pfd = {0};
        pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
        pfd.nVersion = 1;
        pfd.dwFlags = PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW;
        pfd.iPixelType = PFD_TYPE_RGBA;
        pfd.cColorBits = 32;
        pfd.cDepthBits = 32;
        pfd.iLayerType = PFD_MAIN_PLANE;
        const int pixel_format = ChoosePixelFormat(dc, &pfd);
        ASSERT(pixel_format, "ChoosePixelFormat failure");

        const BOOL set_pixel_format_success = SetPixelFormat(dc, pixel_format, &pfd);
        ASSERT(set_pixel_format_success, "SetPixelFormat failure");

        HGLRC temp_context = wglCreateContext(dc);
        ASSERT(temp_context != NULL, "Could not make opengl context.");

        const BOOL wgl_make_current_temp_success = wglMakeCurrent(dc, temp_context);
        ASSERT(wgl_make_current_temp_success == TRUE, "wglMakeCurrent failed.");

        g_opengl_state->wglCreateContextAttribsARB = (fnptr_wglCreateContextAttribsARB)load_gl_fn(opengl32_dll_module, "wglCreateContextAttribsARB");

        int attribs[] =
        {
            WGL_CONTEXT_MAJOR_VERSION_ARB, 4,
            WGL_CONTEXT_MINOR_VERSION_ARB, 4,
            WGL_CONTEXT_FLAGS_ARB, 0,
            0
        };
        g_opengl_state->gl_context = g_opengl_state->wglCreateContextAttribsARB(dc, 0, attribs);
        ASSERT(g_opengl_state->gl_context != NULL, "Could not make opengl context.");

        const BOOL wgl_make_current_null_success = wglMakeCurrent(NULL, NULL);
        ASSERT(wgl_make_current_null_success == TRUE, "wglMakeCurrent failed.");

        const BOOL wgl_make_currentxt_success = wglDeleteContext(temp_context);
        ASSERT(wgl_make_currentxt_success == TRUE, "wglDeleteContext failed.");

        const BOOL wgl_make_current_final_success = wglMakeCurrent(dc, g_opengl_state->gl_context);
        ASSERT(wgl_make_current_final_success == TRUE, "wglMakeCurrent failed.");

        // Load OpenGL functions.
        g_opengl_state->glGetError = (fnptr_glGetError)load_gl_fn(opengl32_dll_module, "glGetError");
        g_opengl_state->glGenVertexArrays = (fnptr_glGenVertexArrays)load_gl_fn(opengl32_dll_module, "glGenVertexArrays");
        g_opengl_state->glBindVertexArray = (fnptr_glBindVertexArray)load_gl_fn(opengl32_dll_module, "glBindVertexArray");
        g_opengl_state->glGenBuffers = (fnptr_glGenBuffers)load_gl_fn(opengl32_dll_module, "glGenBuffers");
        g_opengl_state->glBindBuffer = (fnptr_glBindBuffer)load_gl_fn(opengl32_dll_module, "glBindBuffer");
        g_opengl_state->glBufferData = (fnptr_glBufferData)load_gl_fn(opengl32_dll_module, "glBufferData");
        g_opengl_state->glVertexAttribPointer = (fnptr_glVertexAttribPointer)load_gl_fn(opengl32_dll_module, "glVertexAttribPointer");
        g_opengl_state->glEnableVertexAttribArray = (fnptr_glEnableVertexAttribArray)load_gl_fn(opengl32_dll_module, "glEnableVertexAttribArray");
        g_opengl_state->glCreateShader = (fnptr_glCreateShader)load_gl_fn(opengl32_dll_module, "glCreateShader");
        g_opengl_state->glShaderSource = (fnptr_glShaderSource)load_gl_fn(opengl32_dll_module, "glShaderSource");
        g_opengl_state->glCompileShader = (fnptr_glCompileShader)load_gl_fn(opengl32_dll_module, "glCompileShader");
        g_opengl_state->glGetShaderiv = (fnptr_glGetShaderiv)load_gl_fn(opengl32_dll_module, "glGetShaderiv");
        g_opengl_state->glGetShaderInfoLog = (fnptr_glGetShaderInfoLog)load_gl_fn(opengl32_dll_module, "glGetShaderInfoLog");
        g_opengl_state->glCreateProgram = (fnptr_glCreateProgram)load_gl_fn(opengl32_dll_module, "glCreateProgram");
        g_opengl_state->glAttachShader = (fnptr_glAttachShader)load_gl_fn(opengl32_dll_module, "glAttachShader");
        g_opengl_state->glLinkProgram = (fnptr_glLinkProgram)load_gl_fn(opengl32_dll_module, "glLinkProgram");
        g_opengl_state->glGetProgramiv = (fnptr_glGetProgramiv)load_gl_fn(opengl32_dll_module, "glGetProgramiv");
        g_opengl_state->glGetProgramInfoLog = (fnptr_glGetProgramInfoLog)load_gl_fn(opengl32_dll_module, "glGetProgramInfoLog");
        g_opengl_state->glDeleteShader = (fnptr_glDeleteShader)load_gl_fn(opengl32_dll_module, "glDeleteShader");
        g_opengl_state->glUseProgram = (fnptr_glUseProgram)load_gl_fn(opengl32_dll_module, "glUseProgram");
        g_opengl_state->glBufferSubData = (fnptr_glBufferSubData)load_gl_fn(opengl32_dll_module, "glBufferSubData");
        g_opengl_state->glDrawElementsInstanced = (fnptr_glDrawElementsInstanced)load_gl_fn(opengl32_dll_module, "glDrawElementsInstanced");

        glViewport(0, 0, (GLint)window_width, (GLint)window_height);

        // Vertex array object
        CALL_GL(glGenVertexArrays, 1, &g_opengl_state->vertex_array_object);
        CALL_GL(glBindVertexArray, g_opengl_state->vertex_array_object);

        // Index buffer object
        CALL_GL(glGenBuffers, 1, &g_opengl_state->index_buffer_object);
        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, g_opengl_state->index_buffer_object);
        CALL_GL(glBufferData, GL_ELEMENT_ARRAY_BUFFER, INDEX_ARRAY_BYTES, 0, GL_DYNAMIC_DRAW);

        // Vertex buffer: X
        CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_x);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_x);
        CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
        CALL_GL(glVertexAttribPointer, 0, 1, GL_FLOAT, GL_FALSE, 0, 0);
        CALL_GL(glEnableVertexAttribArray, 0);

        // Vertex buffer: Y
        CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_y);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_y);
        CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
        CALL_GL(glVertexAttribPointer, 1, 1, GL_FLOAT, GL_FALSE, 0, 0);
        CALL_GL(glEnableVertexAttribArray, 1);

        // Vertex buffer: Z
        CALL_GL(glGenBuffers, 1, &g_opengl_state->vertex_buffer_object_z);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_z);
        CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
        CALL_GL(glVertexAttribPointer, 2, 1, GL_FLOAT, GL_FALSE, 0, 0);
        CALL_GL(glEnableVertexAttribArray, 2);

        CALL_GL(glBindVertexArray, 0);
        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

        // Shaders
        const char* vertex_source =
            "#version 440 core\n"
            "layout (location = 0) in float a_x;\n"
            "layout (location = 1) in float a_y;\n"
            "layout (location = 2) in float a_z;\n"
            "void main()\n"
            "{\n"
            "    vec3 pos = vec3(a_x, a_y, a_z);\n"
            "    gl_Position = vec4(pos, 1.0f);\n"
            "}\n"
            "\n";
        const char* fragment_source =
            "#version 440 core\n"
            "out vec4 result_frag_color;\n"
            "void main()\n"
            "{\n"
            "    result_frag_color = vec4(1.0f, 0.0f, 1.0f, 1.0f);\n"
            "}\n"
            "\n";

        GLsizei debug_info_len;
        char debug_info_buf[512];
        GLint shader_compile_success;

        GLuint vertex_shader;
        CALL_GL_RET(&vertex_shader, GLuint, glCreateShader, GL_VERTEX_SHADER);
        CALL_GL(glShaderSource, vertex_shader, 1, &vertex_source, NULL);
        CALL_GL(glCompileShader, vertex_shader);
        CALL_GL(glGetShaderiv, vertex_shader, GL_COMPILE_STATUS, &shader_compile_success);
        if(!shader_compile_success)
        {
            CALL_GL(glGetShaderInfoLog, vertex_shader, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
            WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
            ASSERT(0, "Failed to compile shader.");
        }

        GLuint fragment_shader;
        CALL_GL_RET(&fragment_shader, GLuint, glCreateShader, GL_FRAGMENT_SHADER);
        CALL_GL(glShaderSource, fragment_shader, 1, &fragment_source, NULL);
        CALL_GL(glCompileShader, fragment_shader);
        CALL_GL(glGetShaderiv, fragment_shader, GL_COMPILE_STATUS, &shader_compile_success);
        if(!shader_compile_success)
        {
            CALL_GL(glGetShaderInfoLog, vertex_shader, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
            WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
            ASSERT(0, "Failed to compile shader.");
        }

        CALL_GL_RET(&g_opengl_state->shader_program, GLuint, glCreateProgram);
        CALL_GL(glAttachShader, g_opengl_state->shader_program, vertex_shader);
        CALL_GL(glAttachShader, g_opengl_state->shader_program, fragment_shader);
        CALL_GL(glLinkProgram, g_opengl_state->shader_program);
        
        CALL_GL(glGetProgramiv, g_opengl_state->shader_program, GL_LINK_STATUS, &shader_compile_success);
        if(!shader_compile_success)
        {
            CALL_GL(glGetProgramInfoLog, g_opengl_state->shader_program, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
            WriteFile(g_log_file, debug_info_buf, debug_info_len, NULL, NULL);
            ASSERT(0, "Failed to link shader.");
        }
        
        CALL_GL(glDeleteShader, vertex_shader);
        CALL_GL(glDeleteShader, fragment_shader); 

    }

    // Main loop
    while(1)
    {
        MSG msg;
        while(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
        {
            // https://learn.microsoft.com/en-us/windows/win32/learnwin32/closing-the-window?redirectedfrom=MSDN
            if(msg.message == WM_QUIT)
            {
                ExitProcess(0);
            }
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }

        if(is_key_down(KB_ESCAPE))
        {
            ExitProcess(0);
        }



        float vx[] = {-0.5f, 0.0f, 0.5f};
        float vy[] = {-0.5f, 0.5f, -0.5f};
        float vz[] = {0.0f, 0.0f, -0.0f};
        u32 indices[] = {0, 1, 2};

        glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        CALL_GL(glUseProgram, g_opengl_state->shader_program);

        //{
        //    mat4 m = view_m_world(current_cam);
        //    GLint loc = glGetUniformLocation(G_graphics_state->batch_voxel_shader_program, "m_view");
        //    glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
        //    if(loc == -1) assert(false);
        //}
        //{
        //    mat4 m = clip_m_view(current_cam);
        //    GLint loc = glGetUniformLocation(G_graphics_state->batch_voxel_shader_program, "m_proj");
        //    glUniformMatrix4fv(loc, 1, true, &(m[0][0]));
        //    if(loc == -1) assert(false);
        //}


        CALL_GL(glBindVertexArray, g_opengl_state->vertex_array_object);


        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_x);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(vx), vx);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_y);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(vy), vy);

        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->vertex_buffer_object_z);
        CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(vz), vz);


        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, g_opengl_state->index_buffer_object);
        CALL_GL(glBufferSubData, GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(indices), indices);


        CALL_GL(glDrawElementsInstanced, GL_TRIANGLES, ARRAY_COUNT(indices), GL_UNSIGNED_INT, 0, 1/*batch_size*/);


        const HDC dc = GetDC(g_opengl_state->hwnd);
        BOOL swap_buffers_success = SwapBuffers(dc);
        ASSERT(swap_buffers_success, "SwapBuffers failed.");
    }

    CloseHandle(g_log_file);

    ExitProcess(0);
}

