
#include "platform.h"
#include "graphics.h"
#include "game_state.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <GL/GL.h>
#include "wglext.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Globals

struct Platform_Win32* g_platform;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// OpenGL / WGL definitions

// Definitions copied from glad library.
typedef signed long long int GLsizeiptr;
typedef signed long long int GLintptr;
typedef char GLchar;
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_ARRAY_BUFFER 0x8892
#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#define GL_CLAMP_TO_EDGE 0x812F


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
typedef void (*fnptr_glDrawArrays)(GLenum mode, GLint first, GLsizei count);
typedef void (*fnptr_glDrawArraysInstanced)(GLenum mode, GLint first, GLsizei count, GLsizei instancecount);
typedef void (*fnptr_glDrawElementsInstanced)(GLenum mode, GLsizei count, GLenum type, const void * indices, GLsizei instancecount);
typedef GLint (*fnptr_glGetUniformLocation)(GLuint program, const GLchar *name);
typedef void (*fnptr_glUniformMatrix4fv)(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (*fnptr_glUniformMatrix4fv)(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (*fnptr_glVertexAttribDivisor)(GLuint index, GLuint divisor);
typedef void (*fnptr_glGenTextures)(GLsizei n, GLuint * textures);
typedef void (*fnptr_glBindTexture)(GLenum target, GLuint texture);
typedef void (*fnptr_glTexParameteri)(GLenum target, GLenum pname, GLint param);
typedef void (*fnptr_glTexImage2D)(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const void * data);
typedef void (*fnptr_glTexSubImage2D)(GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum format, GLenum type, const void * pixels);


#define VERTEX_ARRAY_BYTES MB(10)
#define INDEX_ARRAY_BYTES MB(10)

struct OpenGLState
{
    HWND hwnd;
    u32_m screen_width;
    u32_m screen_height;

    HGLRC gl_context;
    GLuint last_gl_error;

    GLuint vertex_array_object;
    GLuint index_buffer_object;
    GLuint vertex_buffer_object_vx;
    GLuint vertex_buffer_object_vy;
    GLuint vertex_buffer_object_vz;
    GLuint vertex_buffer_object_nx;
    GLuint vertex_buffer_object_ny;
    GLuint vertex_buffer_object_nz;
    GLuint instanced_vertex_buffer_object_offset_x;
    GLuint instanced_vertex_buffer_object_offset_y;
    GLuint instanced_vertex_buffer_object_offset_z;
    GLuint shader_program;

    GLuint debug_line_vertex_array_object;
    GLuint debug_line_vertex_buffer_object_vertices;
    GLuint debug_line_instanced_vertex_buffer_object_start_x;
    GLuint debug_line_instanced_vertex_buffer_object_start_y;
    GLuint debug_line_instanced_vertex_buffer_object_start_z;
    GLuint debug_line_instanced_vertex_buffer_object_end_x;
    GLuint debug_line_instanced_vertex_buffer_object_end_y;
    GLuint debug_line_instanced_vertex_buffer_object_end_z;
    GLuint debug_line_instanced_vertex_buffer_object_color_r;
    GLuint debug_line_instanced_vertex_buffer_object_color_g;
    GLuint debug_line_instanced_vertex_buffer_object_color_b;
    GLuint debug_line_shader_program;

    GLuint debug_sphere_vertex_array_object;
    GLuint debug_sphere_vertex_buffer_object_vertices;
    GLuint debug_sphere_instanced_vertex_buffer_object_pos_x;
    GLuint debug_sphere_instanced_vertex_buffer_object_pos_y;
    GLuint debug_sphere_instanced_vertex_buffer_object_pos_z;
    GLuint debug_sphere_instanced_vertex_buffer_object_radius;
    GLuint debug_sphere_instanced_vertex_buffer_object_color_r;
    GLuint debug_sphere_instanced_vertex_buffer_object_color_g;
    GLuint debug_sphere_instanced_vertex_buffer_object_color_b;
    GLuint debug_sphere_shader_program;

    GLuint fullscreen_quad_vertex_array_object;
    GLuint fullscreen_quad_vertex_buffer_object_pos;
    GLuint fullscreen_quad_texture;
    GLuint textured_quad_shader_program;

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
    fnptr_glDrawArrays glDrawArrays;
    fnptr_glDrawArraysInstanced glDrawArraysInstanced;
    fnptr_glDrawElementsInstanced glDrawElementsInstanced;
    fnptr_glGetUniformLocation glGetUniformLocation;
    fnptr_glUniformMatrix4fv glUniformMatrix4fv;
    fnptr_glVertexAttribDivisor glVertexAttribDivisor;
    fnptr_glGenTextures glGenTextures;
    fnptr_glBindTexture glBindTexture;
    fnptr_glTexParameteri glTexParameteri;
    fnptr_glTexImage2D glTexImage2D;
    fnptr_glTexSubImage2D glTexSubImage2D;
};

INTERNAL f32 s_cube_mesh_vx[] = { -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f };
INTERNAL f32 s_cube_mesh_vy[] = { -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f };
INTERNAL f32 s_cube_mesh_vz[] = { -0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f };
INTERNAL u32 s_cube_mesh_indices[] = 
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Input definitions

#define MAX_MOUSE_KEYS 32
struct InputState
{
    u32_m key[NUM_KEYBOARD_KEYS];
    u32_m mouse_key[MAX_MOUSE_KEYS];

    u32_m last_key[NUM_KEYBOARD_KEYS];
    u32_m last_mouse_key[MAX_MOUSE_KEYS];

    u32_m fps_mode; 
    s32_m mouse_screen_pos_x;
    s32_m mouse_screen_pos_y;
    s32_m mouse_screen_dx;
    s32_m mouse_screen_dy;
};

// Maps Win32 VK code to KeyboardKey enum.
INTERNAL u32 win32_vk_code_to_keyboard_key[] = {
    // https://learn.microsoft.com/en-us/windows/win32/inputdev/virtual-key-codes
    /* -              0x00 */  KB_NOT_SUPPORTED,
    /* VK_LBUTTON     0x01 */  KB_NOT_SUPPORTED,
    /* VK_RBUTTON     0x02 */  KB_NOT_SUPPORTED,
    /* VK_CANCEL      0x03 */  KB_NOT_SUPPORTED,
    /* VK_MBUTTON     0x04 */  KB_NOT_SUPPORTED,
    /* VK_XBUTTON1    0x05 */  KB_NOT_SUPPORTED,
    /* VK_XBUTTON2    0x06 */  KB_NOT_SUPPORTED,
    /* -              0x07 */  KB_NOT_SUPPORTED,
    /* VK_BACK        0x08 */  KB_NOT_SUPPORTED,
    /* VK_TAB         0x09 */  KB_NOT_SUPPORTED,
    /* -              0x0A */  KB_NOT_SUPPORTED,
    /* -              0x0B */  KB_NOT_SUPPORTED,
    /* VK_CLEAR       0x0C */  KB_NOT_SUPPORTED,
    /* VK_RETURN      0x0D */  KB_NOT_SUPPORTED,
    /* -              0x0E */  KB_NOT_SUPPORTED,
    /* -              0x0F */  KB_NOT_SUPPORTED,

    /* VK_SHIFT       0x10 */  KB_NOT_SUPPORTED,
    /* VK_CONTROL     0x11 */  KB_LCTRL,
    /* VK_MENU        0x12 */  KB_NOT_SUPPORTED,
    /* VK_PAUSE       0x13 */  KB_NOT_SUPPORTED,
    /* VK_CAPITAL     0x14 */  KB_NOT_SUPPORTED,
    /* VK_KANA        0x15 */  KB_NOT_SUPPORTED,
    /* VK_IME_ON      0x16 */  KB_NOT_SUPPORTED,
    /* VK_JUNJA       0x17 */  KB_NOT_SUPPORTED,
    /* VK_FINAL       0x18 */  KB_NOT_SUPPORTED,
    /* VK_HANJA       0x19 */  KB_NOT_SUPPORTED,
    /* VK_IME_OFF     0x1A */  KB_NOT_SUPPORTED,
    /* VK_ESCAPE      0x1B */  KB_ESCAPE,
    /* VK_CONVERT     0x1C */  KB_NOT_SUPPORTED,
    /* VK_NONCONVERT  0x1D */  KB_NOT_SUPPORTED,
    /* VK_ACCEPT      0x1E */  KB_NOT_SUPPORTED,
    /* VK_MODECHANGE  0x1F */  KB_NOT_SUPPORTED,

    /* VK_SPACE       0x20 */  KB_SPACE,
    /* VK_PRIOR       0x21 */  KB_NOT_SUPPORTED,
    /* VK_NEXT        0x22 */  KB_NOT_SUPPORTED,
    /* VK_END         0x23 */  KB_NOT_SUPPORTED,
    /* VK_HOME        0x24 */  KB_NOT_SUPPORTED,
    /* VK_LEFT        0x25 */  KB_NOT_SUPPORTED,
    /* VK_UP          0x26 */  KB_NOT_SUPPORTED,
    /* VK_RIGHT       0x27 */  KB_NOT_SUPPORTED,
    /* VK_DOWN        0x28 */  KB_NOT_SUPPORTED,
    /* VK_SELECT      0x29 */  KB_NOT_SUPPORTED,
    /* VK_PRINT       0x2A */  KB_NOT_SUPPORTED,
    /* VK_EXECUTE     0x2B */  KB_NOT_SUPPORTED,
    /* VK_SNAPSHOT    0x2C */  KB_NOT_SUPPORTED,
    /* VK_INSERT      0x2D */  KB_NOT_SUPPORTED,
    /* VK_DELETE      0x2E */  KB_NOT_SUPPORTED,
    /* VK_HELP        0x2F */  KB_NOT_SUPPORTED,

    /* 0 key          0x30 */  KB_0,
    /* 1 key          0x31 */  KB_1,
    /* 2 key          0x32 */  KB_2,
    /* 3 key          0x33 */  KB_3,
    /* 4 key          0x34 */  KB_4,
    /* 5 key          0x35 */  KB_5,
    /* 6 key          0x36 */  KB_6,
    /* 7 key          0x37 */  KB_7,
    /* 8 key          0x38 */  KB_8,
    /* 9 key          0x39 */  KB_9,
    /* Undefined      0x3A */  KB_NOT_SUPPORTED,
    /* Undefined      0x3B */  KB_NOT_SUPPORTED,
    /* Undefined      0x3C */  KB_NOT_SUPPORTED,
    /* Undefined      0x3D */  KB_NOT_SUPPORTED,
    /* Undefined      0x3E */  KB_NOT_SUPPORTED,
    /* Undefined      0x3F */  KB_NOT_SUPPORTED,

    /* Undefined      0x40 */  KB_NOT_SUPPORTED,
    /* A key          0x41 */  KB_A,
    /* B key          0x42 */  KB_B,
    /* C key          0x43 */  KB_C,
    /* D key          0x44 */  KB_D,
    /* E key          0x45 */  KB_E,
    /* F key          0x46 */  KB_F,
    /* G key          0x47 */  KB_G,
    /* H key          0x48 */  KB_H,
    /* I key          0x49 */  KB_I,
    /* J key          0x4A */  KB_J,
    /* K key          0x4B */  KB_K,
    /* L key          0x4C */  KB_L,
    /* M key          0x4D */  KB_M,
    /* N key          0x4E */  KB_N,
    /* O key          0x4F */  KB_O,

    /* P key          0x50 */  KB_P,
    /* Q key          0x51 */  KB_Q,
    /* R key          0x52 */  KB_R,
    /* S key          0x53 */  KB_S,
    /* T key          0x54 */  KB_T,
    /* U key          0x55 */  KB_U,
    /* V key          0x56 */  KB_V,
    /* W key          0x57 */  KB_W,
    /* X key          0x58 */  KB_X,
    /* Y key          0x59 */  KB_Y,
    /* Z key          0x5A */  KB_Z,
    /* VK_LWIN        0x5B */  KB_NOT_SUPPORTED,
    /* VK_RWIN        0x5C */  KB_NOT_SUPPORTED,
    /* VK_APPS        0x5D */  KB_NOT_SUPPORTED,
    /* -              0x5E */  KB_NOT_SUPPORTED,
    /* VK_SLEEP       0x5F */  KB_NOT_SUPPORTED,

    /* VK_NUMPAD0     0x60 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD1     0x61 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD2     0x62 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD3     0x63 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD4     0x64 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD5     0x65 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD6     0x66 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD7     0x67 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD8     0x68 */  KB_NOT_SUPPORTED,
    /* VK_NUMPAD9     0x69 */  KB_NOT_SUPPORTED,
    /* VK_MULTIPLY    0x6A */  KB_NOT_SUPPORTED,
    /* VK_ADD         0x6B */  KB_NOT_SUPPORTED,
    /* VK_SEPARATOR   0x6C */  KB_NOT_SUPPORTED,
    /* VK_SUBTRACT    0x6D */  KB_NOT_SUPPORTED,
    /* VK_DECIMAL     0x6E */  KB_NOT_SUPPORTED,
    /* VK_DIVIDE      0x6F */  KB_NOT_SUPPORTED,

    /* VK_F1          0x70 */  KB_NOT_SUPPORTED,
    /* VK_F2          0x71 */  KB_NOT_SUPPORTED,
    /* VK_F3          0x72 */  KB_NOT_SUPPORTED,
    /* VK_F4          0x73 */  KB_NOT_SUPPORTED,
    /* VK_F5          0x74 */  KB_NOT_SUPPORTED,
    /* VK_F6          0x75 */  KB_NOT_SUPPORTED,
    /* VK_F7          0x76 */  KB_NOT_SUPPORTED,
    /* VK_F8          0x77 */  KB_NOT_SUPPORTED,
    /* VK_F9          0x78 */  KB_NOT_SUPPORTED,
    /* VK_F10         0x79 */  KB_NOT_SUPPORTED,
    /* VK_F11         0x7A */  KB_NOT_SUPPORTED,
    /* VK_F12         0x7B */  KB_NOT_SUPPORTED,
    /* VK_F13         0x7C */  KB_NOT_SUPPORTED,
    /* VK_F14         0x7D */  KB_NOT_SUPPORTED,
    /* VK_F15         0x7E */  KB_NOT_SUPPORTED,
    /* VK_F16         0x7F */  KB_NOT_SUPPORTED,

    /* VK_F17         0x80 */  KB_NOT_SUPPORTED,
    /* VK_F18         0x81 */  KB_NOT_SUPPORTED,
    /* VK_F19         0x82 */  KB_NOT_SUPPORTED,
    /* VK_F20         0x83 */  KB_NOT_SUPPORTED,
    /* VK_F21         0x84 */  KB_NOT_SUPPORTED,
    /* VK_F22         0x85 */  KB_NOT_SUPPORTED,
    /* VK_F23         0x86 */  KB_NOT_SUPPORTED,
    /* VK_F24         0x87 */  KB_NOT_SUPPORTED,
    /* -              0x88 */  KB_NOT_SUPPORTED,
    /* -              0x89 */  KB_NOT_SUPPORTED,
    /* -              0x8A */  KB_NOT_SUPPORTED,
    /* -              0x8B */  KB_NOT_SUPPORTED,
    /* -              0x8C */  KB_NOT_SUPPORTED,
    /* -              0x8D */  KB_NOT_SUPPORTED,
    /* -              0x8E */  KB_NOT_SUPPORTED,
    /* -              0x8F */  KB_NOT_SUPPORTED,

    /* VK_NUMLOCK     0x90 */  KB_NOT_SUPPORTED,
    /* VK_SCROLL      0x91 */  KB_NOT_SUPPORTED,
    /* -              0x92 */  KB_NOT_SUPPORTED,
    /* -              0x93 */  KB_NOT_SUPPORTED,
    /* -              0x94 */  KB_NOT_SUPPORTED,
    /* -              0x95 */  KB_NOT_SUPPORTED,
    /* -              0x96 */  KB_NOT_SUPPORTED,
    /* -              0x97 */  KB_NOT_SUPPORTED,
    /* -              0x98 */  KB_NOT_SUPPORTED,
    /* -              0x99 */  KB_NOT_SUPPORTED,
    /* -              0x9A */  KB_NOT_SUPPORTED,
    /* -              0x9B */  KB_NOT_SUPPORTED,
    /* -              0x9C */  KB_NOT_SUPPORTED,
    /* -              0x9D */  KB_NOT_SUPPORTED,
    /* -              0x9E */  KB_NOT_SUPPORTED,
    /* -              0x9F */  KB_NOT_SUPPORTED,

    /* VK_LSHIFT      0xA0 */  KB_NOT_SUPPORTED,
    /* VK_RSHIFT      0xA1 */  KB_NOT_SUPPORTED,
    /* VK_LCONTROL    0xA2 */  KB_LCTRL,
    /* VK_RCONTROL    0xA3 */  KB_NOT_SUPPORTED,
};
_Static_assert(ARRAY_COUNT(win32_vk_code_to_keyboard_key) == 0xA4, "VK map is not the expected size.");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Platform Win32 definitions

struct Platform_Win32
{
    HANDLE log_file;
    HWND hwnd;
    u64_m clock_freq;

    struct OpenGLState opengl_state;
    struct InputState input_state;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Platform Win32 impl
void assert_fn(const char* file, int line, u64 c, const char* msg)
{
    if(!c)
    {
        const StringBuf32 line_str = u32_to_StringBuf32(line);
        const char* line_cstr = StringBuf32_to_cstr(&line_str);
        WriteFile(g_platform->log_file, file, (int)cstr_len(file), NULL, NULL);
        WriteFile(g_platform->log_file, ":", 1, NULL, NULL);
        WriteFile(g_platform->log_file, line_cstr, (int)cstr_len(line_cstr), NULL, NULL);
        WriteFile(g_platform->log_file, "  ", 1, NULL, NULL);
        WriteFile(g_platform->log_file, msg, (int)cstr_len(msg), NULL, NULL);
        FlushFileBuffers(g_platform->log_file);
        DebugBreak();
    }
}

struct MemoryArena
{
    u8* const end;
    u8* const start;
    u8_m* cur;
};

struct MemoryArena memory_arena_init(void* base, u64 cap)
{
    ASSERT(base, "Memory arena init with NULL base.");
    struct MemoryArena result = {
        .start = base,
        .end = (u8*)base + cap,
        .cur = (u8_m*)base,
    };
    return result;
}

void memory_arena_reset(struct MemoryArena* arena)
{
    arena->cur = (u8_m*)arena->start;
}

void* memory_arena_allocate(const char* file, s32 line, struct MemoryArena* arena, u64 req_bytes)
{
    (void)file;
    (void)line;
    // Align to 64 byte boundary.
    u64 bytes = (req_bytes + 63ULL) & ~63ULL;
    ASSERT(arena->cur + bytes <= arena->end, "Memory arena allocation overflow.");
    void* result = arena->cur;
    arena->cur += bytes;
    return result;
}

void* memory_arena_allocate_zeroed(const char* file, s32 line, struct MemoryArena* arena, u64 req_bytes)
{
    (void)file;
    (void)line;
    void* result = memory_arena_allocate(file, line, arena, req_bytes);
    memset(result, 0, req_bytes);
    return result;
}

LRESULT CALLBACK WindowProc(HWND window, UINT message, WPARAM wParam, LPARAM lParam)
{
    LRESULT result = 0; 

    struct Platform_Win32* platform = g_platform;

    switch (message)
    {
        case WM_KEYDOWN: 
        {
            u32 vk = (u32)wParam;
            ASSERT(vk < ARRAY_COUNT(win32_vk_code_to_keyboard_key), "Invalid VK code");
            u32 k = win32_vk_code_to_keyboard_key[vk];
            platform->input_state.key[k] = 1;
            break;
        }

        case WM_KEYUP:
        {
            u32 vk = (u32)wParam;
            ASSERT(vk < ARRAY_COUNT(win32_vk_code_to_keyboard_key), "Invalid VK code");
            u32 k = win32_vk_code_to_keyboard_key[vk];
            platform->input_state.key[k] = 0;
            break;
        }

        case WM_LBUTTONDOWN:
        {
            platform->input_state.mouse_key[0] = 1;
            break;
        }

        case WM_LBUTTONUP:
        {
            platform->input_state.mouse_key[0] = 0;
            break;
        }

        case WM_RBUTTONDOWN:
        {
            platform->input_state.mouse_key[1] = 1;
            break;
        }

        case WM_RBUTTONUP:
        {
            platform->input_state.mouse_key[1] = 0;
            break;
        }

        case WM_MOUSEMOVE: 
        {
            struct InputState* input_state = &platform->input_state;
            if(input_state->fps_mode)
            {
                RECT clip_rect;
                const BOOL get_client_rect_success = GetWindowRect(platform->hwnd, &clip_rect);
                ASSERT(get_client_rect_success, "GetWindowRect failed.");
                const BOOL clip_cursor_success = ClipCursor(&clip_rect);
                ASSERT(clip_cursor_success, "ClipCursor failed.");
            }
            else
            {
                const BOOL clip_cursor_success = ClipCursor(NULL);
                ASSERT(clip_cursor_success, "ClipCursor failed.");
            }
            break;
        }

        case WM_SIZE:
        {
            const UINT width = LOWORD(lParam);
            const UINT height = HIWORD(lParam);
            platform->opengl_state.screen_width = width;
            platform->opengl_state.screen_height = height;
            glViewport(0, 0, width, height);
            break;
        }

        case WM_KILLFOCUS:
        {
            struct InputState* input_state = &platform->input_state;
            input_state->fps_mode = 0;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// OpenGL / WGL impl

// https://www.khronos.org/opengl/wiki/Load_OpenGL_Functions
INTERNAL void *load_gl_fn(HMODULE opengl32_dll_module, const char *name)
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

#define CALL_GL(fn, ...) \
do { \
    struct OpenGLState* CALL_GL_opengl_state = &g_platform->opengl_state; \
    CALL_GL_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = CALL_GL_opengl_state->glGetError(); \
    CALL_GL_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
} while(0)

#define CALL_GL_RET(out, type, fn, ...) \
do { \
    struct OpenGLState* CALL_GL_opengl_state = &g_platform->opengl_state; \
    type MACRO_r = CALL_GL_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = CALL_GL_opengl_state->glGetError(); \
    CALL_GL_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
    *(out) = MACRO_r; \
} while(0)


INTERNAL GLuint make_shader_program(const char* vertex_source, const char* fragment_source)
{
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
        WriteFile(g_platform->log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to compile vertex shader.");
    }

    GLuint fragment_shader;
    CALL_GL_RET(&fragment_shader, GLuint, glCreateShader, GL_FRAGMENT_SHADER);
    CALL_GL(glShaderSource, fragment_shader, 1, &fragment_source, NULL);
    CALL_GL(glCompileShader, fragment_shader);
    CALL_GL(glGetShaderiv, fragment_shader, GL_COMPILE_STATUS, &shader_compile_success);
    if(!shader_compile_success)
    {
        CALL_GL(glGetShaderInfoLog, fragment_shader, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
        WriteFile(g_platform->log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to compile fragment shader.");
    }

    u32_m result_program;
    CALL_GL_RET(&result_program, GLuint, glCreateProgram);
    CALL_GL(glAttachShader, result_program, vertex_shader);
    CALL_GL(glAttachShader, result_program, fragment_shader);
    CALL_GL(glLinkProgram, result_program);
    
    CALL_GL(glGetProgramiv, result_program, GL_LINK_STATUS, &shader_compile_success);
    if(!shader_compile_success)
    {
        CALL_GL(glGetProgramInfoLog, result_program, sizeof(debug_info_buf), &debug_info_len, debug_info_buf);
        WriteFile(g_platform->log_file, debug_info_buf, debug_info_len, NULL, NULL);
        ASSERT(0, "Failed to link shader.");
    }
    
    CALL_GL(glDeleteShader, vertex_shader);
    CALL_GL(glDeleteShader, fragment_shader); 

    return result_program;
}

f32 get_screen_aspect_ratio()
{
    const struct OpenGLState* opengl_state = &g_platform->opengl_state;
    return (float)opengl_state->screen_height / (float)opengl_state->screen_width;
}

void draw_Mesh_64(const struct Mesh_64* mesh, const mtx4x4* camera_and_clip_mtx, u64 num, f32* pos_x, f32* pos_y, f32* pos_z)
{
    struct OpenGLState* opengl_state = &g_platform->opengl_state;
    CALL_GL(glUseProgram, opengl_state->shader_program);

    GLint loc;
    CALL_GL_RET(&loc, GLint, glGetUniformLocation, opengl_state->shader_program, "m_mvp");
    CALL_GL(glUniformMatrix4fv, loc, 1, 1, &camera_and_clip_mtx->m[0]);
    ASSERT(loc != -1, "Failed to bind uniform.");

    CALL_GL(glBindVertexArray, opengl_state->vertex_array_object);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vx);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vx);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vy);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vy);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vz);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vz);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nx);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->nx);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_ny);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->ny);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nz);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->nz);

    ASSERT(num * sizeof(f32) <= VERTEX_ARRAY_BYTES, "Instanced buffer overflow.");
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_x);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_x);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_y);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_y);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_z);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_z);

    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, opengl_state->index_buffer_object);
    CALL_GL(glBufferSubData, GL_ELEMENT_ARRAY_BUFFER, 0, mesh->num_indices * sizeof(u32), mesh->indices);

    CALL_GL(glDrawElementsInstanced, GL_TRIANGLES, mesh->num_indices, GL_UNSIGNED_INT, 0, (GLsizei)num);

    CALL_GL(glBindVertexArray, 0);
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
    CALL_GL(glUseProgram, 0);

}

void draw_Mesh_1M(const struct Mesh_1M* mesh, const mtx4x4* camera_and_clip_mtx, u64 num, f32* pos_x, f32* pos_y, f32* pos_z)
{
    struct OpenGLState* opengl_state = &g_platform->opengl_state;
    CALL_GL(glUseProgram, opengl_state->shader_program);

    GLint loc;
    CALL_GL_RET(&loc, GLint, glGetUniformLocation, opengl_state->shader_program, "m_mvp");
    CALL_GL(glUniformMatrix4fv, loc, 1, 1, &camera_and_clip_mtx->m[0]);
    ASSERT(loc != -1, "Failed to bind uniform.");

    CALL_GL(glBindVertexArray, opengl_state->vertex_array_object);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vx);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vx);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vy);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vy);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vz);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->vz);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nx);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->nx);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_ny);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->ny);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nz);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, mesh->num_vertices * sizeof(f32), mesh->nz);

    ASSERT(num * sizeof(f32) <= VERTEX_ARRAY_BYTES, "Instanced buffer overflow.");
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_x);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_x);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_y);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_y);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_z);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num * sizeof(f32), pos_z);

    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, opengl_state->index_buffer_object);
    CALL_GL(glBufferSubData, GL_ELEMENT_ARRAY_BUFFER, 0, mesh->num_indices * sizeof(u32), mesh->indices);

    CALL_GL(glDrawElementsInstanced, GL_TRIANGLES, mesh->num_indices, GL_UNSIGNED_INT, 0, (GLsizei)num);

    CALL_GL(glBindVertexArray, 0);
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
    CALL_GL(glUseProgram, 0);

}

void debug_draw_line(const mtx4x4* camera_and_clip_mtx, v3 a, v3 b, v3 c)
{
    struct OpenGLState* opengl_state = &g_platform->opengl_state;
    CALL_GL(glUseProgram, opengl_state->debug_line_shader_program);
    
    GLint debug_line_loc;
    CALL_GL_RET(&debug_line_loc, GLint, glGetUniformLocation, opengl_state->shader_program, "m_mvp");
    CALL_GL(glUniformMatrix4fv, debug_line_loc, 1, 1, &camera_and_clip_mtx->m[0]);
    ASSERT(debug_line_loc != -1, "Failed to bind uniform.");
    
    CALL_GL(glBindVertexArray, opengl_state->debug_line_vertex_array_object);

    f32 line_vertices[6] = {
        0.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f,
    };
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_vertex_buffer_object_vertices);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(line_vertices), line_vertices);

    u32 num_lines = 1;
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), a.m + 0);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), a.m + 1);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), a.m + 2);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), b.m + 0);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), b.m + 1);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), b.m + 2);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), c.m + 0);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), c.m + 1);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_lines * sizeof(float), c.m + 2);
    
    CALL_GL(glDrawArraysInstanced, GL_LINES, 0, 2, num_lines);

    CALL_GL(glBindVertexArray, 0);
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
    CALL_GL(glUseProgram, 0);
}

void debug_draw_sphere(const mtx4x4* proj_mtx, const mtx4x4* cam_mtx, v3 pos, f32 r, v3 c)
{
    struct OpenGLState* opengl_state = &g_platform->opengl_state;
    CALL_GL(glUseProgram, opengl_state->debug_sphere_shader_program);
    
    GLint loc;
    CALL_GL_RET(&loc, GLint, glGetUniformLocation, opengl_state->debug_sphere_shader_program, "m_cam");
    CALL_GL(glUniformMatrix4fv, loc, 1, 1, &cam_mtx->m[0]);
    ASSERT(loc != -1, "Failed to bind uniform.");

    CALL_GL_RET(&loc, GLint, glGetUniformLocation, opengl_state->debug_sphere_shader_program, "m_proj");
    CALL_GL(glUniformMatrix4fv, loc, 1, 1, &proj_mtx->m[0]);
    ASSERT(loc != -1, "Failed to bind uniform.");
    
    CALL_GL(glBindVertexArray, opengl_state->debug_sphere_vertex_array_object);

    /*
      3----4
      |    |
      1----2
    */
    f32 verts[3 * 4] = {
        // X      Y     Z
        -1.0f, -1.0f, 0.0f,
         1.0f, -1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f, 0.0f,
    };
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_vertex_buffer_object_vertices);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(verts), verts);

    u32 num_spheres = 1;
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_x);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), pos.m + 0);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_y);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), pos.m + 1);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_z);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), pos.m + 2);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_radius);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), &r);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_r);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), c.m + 0);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_g);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), c.m + 1);

    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_b);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, num_spheres * sizeof(float), c.m + 2);
    
    CALL_GL(glDrawArraysInstanced, GL_TRIANGLE_STRIP, 0, 4, num_spheres);

    CALL_GL(glBindVertexArray, 0);
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
    CALL_GL(glUseProgram, 0);
}

#if 0
void draw_frame_buffer(struct FrameBuffer* fb)
{
    CALL_GL(glUseProgram, g_opengl_state->textured_quad_shader_program);
    CALL_GL(glTexSubImage2D,
            GL_TEXTURE_2D, // GLenum target
            0, // GLint level
            0, // GLint xoffset
            0, // GLint yoffset
            g_draw_data->frame_buffer_width,  // GLsizei width
            g_draw_data->frame_buffer_height, // GLsizei height
            GL_RGBA,          // GLenum format
            GL_UNSIGNED_BYTE, // GLenum type
            g_draw_data->frame_buffer);
    const GLuint asdferr = g_opengl_state->glGetError();
    ASSERT(asdferr == 0, "failed"); \

    CALL_GL(glBindVertexArray, g_opengl_state->fullscreen_quad_vertex_array_object);

    f32 quad_vertices[] = {
        -1.0f, -1.0f,
         1.0f, -1.0f,
         1.0f,  1.0f,

        -1.0f, -1.0f,
         1.0f,  1.0f,
        -1.0f,  1.0f,
    };
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, g_opengl_state->fullscreen_quad_vertex_buffer_object_pos);
    CALL_GL(glBufferSubData, GL_ARRAY_BUFFER, 0, sizeof(quad_vertices), quad_vertices);

    CALL_GL(glDrawArrays, GL_TRIANGLES, 0, 6);

    CALL_GL(glBindVertexArray, 0);
    CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
    CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
    CALL_GL(glUseProgram, 0);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Time

INTERNAL s64 get_timestamp_us()
{
    // https://learn.microsoft.com/en-us/windows/win32/sysinfo/acquiring-high-resolution-time-stamps
    LARGE_INTEGER cy;
    QueryPerformanceCounter(&cy);
    cy.QuadPart *= 1000000LL;
    cy.QuadPart /= g_platform->clock_freq;
    s64 result_us = cy.QuadPart;
    ASSERT(result_us >= 0, "get_timestamp_us overflow.");
    return result_us;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Input impl


u32 is_key_down(enum KeyboardKey k)
{
    ASSERT(k >= 0, "is_key_down key code underflow");
    ASSERT(k < NUM_KEYBOARD_KEYS, "is_key_down key code overflow");
    struct InputState* input_state = &g_platform->input_state;
    return input_state->key[(u32)k];
}

u32 is_key_toggled_down(enum KeyboardKey k)
{
    ASSERT(k >= 0, "is_key_down key code underflow");
    ASSERT(k < NUM_KEYBOARD_KEYS, "is_key_down key code overflow");
    struct InputState* input_state = &g_platform->input_state;
    return input_state->key[(u32)k] && !input_state->last_key[(u32)k];
}

u32 is_fps_mode()
{
    struct InputState* input_state = &g_platform->input_state;
    return input_state->fps_mode;
}

void get_mouse_delta(s32_m* x, s32_m* y)
{
    struct InputState* input_state = &g_platform->input_state;
    *x = input_state->mouse_screen_dx;
    *y = input_state->mouse_screen_dy;
}

void show_cursor(u32_m show)
{
    // https://learn.microsoft.com/en-us/windows/win32/api/winuser/ns-winuser-cursorinfo
    CURSORINFO ci = { .cbSize = sizeof(CURSORINFO) };
    BOOL get_cursor_info_success = GetCursorInfo(&ci);
    ASSERT(get_cursor_info_success, "GetCursorInfo failed.");
    if(show && ci.flags == 0)
    {
        // https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-showcursor
        s32_m show_cursor_display_counter = ShowCursor(1);
        while(show_cursor_display_counter < 0)
        {
            show_cursor_display_counter = ShowCursor(1);
        }
    }
    else if(!show && ci.flags == 0x1)
    {
        // https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-showcursor
        s32_m show_cursor_display_counter = ShowCursor(0);
        while(show_cursor_display_counter >= 0)
        {
            show_cursor_display_counter = ShowCursor(0);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Main

void do_one_frame(struct GameState* prev_game_state, struct GameState* next_game_state, struct MemoryArena*);
void draw_game_state(struct GameState* game_state, struct MemoryArena* memory_arena);

void WinMainCRTStartup()
{
    // Init application memory.
    void* main_memory_arena_starting_address = (void*)0x100000;
    u64 main_memory_arena_cap = GB(1);
    void* main_memory_arena_storage = VirtualAlloc(
        main_memory_arena_starting_address,
        main_memory_arena_cap,
        MEM_RESERVE | MEM_COMMIT,
        PAGE_READWRITE
    );
    ASSERT(main_memory_arena_storage, "Could not allocate main memory arena.");
    struct MemoryArena main_memory_arena = memory_arena_init(main_memory_arena_storage, main_memory_arena_cap);


    g_platform = (struct Platform_Win32*)MEMORY_ARENA_ALLOCATE_ZEROED(&main_memory_arena, sizeof(struct Platform_Win32));
    struct Platform_Win32* platform = g_platform;


    // Init logging file.
    g_platform->log_file = CreateFileA("log.txt", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    ASSERT(g_platform->log_file != INVALID_HANDLE_VALUE, "Couldn't open log file for writing.");


    // Init engine graphics state.
    {
        struct OpenGLState* opengl_state = &g_platform->opengl_state;
        const HMODULE hinstance = GetModuleHandle(0);

        WNDCLASS window_class = {0};
        window_class.style = CS_HREDRAW|CS_VREDRAW|CS_OWNDC;
        window_class.lpfnWndProc = WindowProc;
        window_class.hInstance = hinstance;
        window_class.hCursor = LoadCursor(NULL, IDC_ARROW);
        window_class.lpszClassName = "Windows Program Class";
        ATOM register_class_result = RegisterClass(&window_class);
        ASSERT(register_class_result, "Register class failed");

        u64 window_width = 1920;
        u64 window_height = 1080;
        u64 window_start_x = (GetSystemMetrics(SM_CXSCREEN) - window_width) / 2;
        u64 window_start_y = (GetSystemMetrics(SM_CYSCREEN) - window_height) / 2;
        platform->hwnd = CreateWindowEx(0,          // Extended style
                window_class.lpszClassName,               // Class name
                "",                                       // Window name
                //WS_OVERLAPPEDWINDOW | WS_VISIBLE,       // Style of the window
                WS_POPUP | WS_VISIBLE,                    // Style of the window
                (int)window_start_x,                      // Initial X position
                (int)window_start_y,                      // Initial Y position
                (int)window_width,                        // Initial width
                (int)window_height,                       // Initial height
                0,                                        // Handle to the window parent
                0,                                        // Handle to a menu
                hinstance,                                // Handle to an instance
                0);
        ASSERT(platform->hwnd != NULL, "Failed to create a window");

        const HDC dc = GetDC(platform->hwnd);
        ASSERT(dc != NULL, "GetDC failed.");

        RECT client_rect;
        const BOOL get_client_rect_success = GetClientRect(platform->hwnd, &client_rect);
        ASSERT(get_client_rect_success, "GetClientRect failed.");
        opengl_state->screen_width = client_rect.right;
        opengl_state->screen_height = client_rect.bottom;

        BOOL clip_cursor_success = ClipCursor(&client_rect);
        ASSERT(clip_cursor_success, "ClipCursor failed.");

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

        opengl_state->wglCreateContextAttribsARB = (fnptr_wglCreateContextAttribsARB)load_gl_fn(opengl32_dll_module, "wglCreateContextAttribsARB");

        int attribs[] =
        {
            WGL_CONTEXT_MAJOR_VERSION_ARB, 3,
            WGL_CONTEXT_MINOR_VERSION_ARB, 3,
            WGL_CONTEXT_FLAGS_ARB, 0,
            0
        };
        opengl_state->gl_context = opengl_state->wglCreateContextAttribsARB(dc, 0, attribs);
        ASSERT(opengl_state->gl_context != NULL, "Could not make opengl context.");

        const BOOL wgl_make_current_null_success = wglMakeCurrent(NULL, NULL);
        ASSERT(wgl_make_current_null_success == TRUE, "wglMakeCurrent failed.");

        const BOOL wgl_make_currentxt_success = wglDeleteContext(temp_context);
        ASSERT(wgl_make_currentxt_success == TRUE, "wglDeleteContext failed.");

        const BOOL wgl_make_current_final_success = wglMakeCurrent(dc, opengl_state->gl_context);
        ASSERT(wgl_make_current_final_success == TRUE, "wglMakeCurrent failed.");

        // Load OpenGL functions.
        opengl_state->glGetError = (fnptr_glGetError)load_gl_fn(opengl32_dll_module, "glGetError");
        opengl_state->glGenVertexArrays = (fnptr_glGenVertexArrays)load_gl_fn(opengl32_dll_module, "glGenVertexArrays");
        opengl_state->glBindVertexArray = (fnptr_glBindVertexArray)load_gl_fn(opengl32_dll_module, "glBindVertexArray");
        opengl_state->glGenBuffers = (fnptr_glGenBuffers)load_gl_fn(opengl32_dll_module, "glGenBuffers");
        opengl_state->glBindBuffer = (fnptr_glBindBuffer)load_gl_fn(opengl32_dll_module, "glBindBuffer");
        opengl_state->glBufferData = (fnptr_glBufferData)load_gl_fn(opengl32_dll_module, "glBufferData");
        opengl_state->glVertexAttribPointer = (fnptr_glVertexAttribPointer)load_gl_fn(opengl32_dll_module, "glVertexAttribPointer");
        opengl_state->glEnableVertexAttribArray = (fnptr_glEnableVertexAttribArray)load_gl_fn(opengl32_dll_module, "glEnableVertexAttribArray");
        opengl_state->glCreateShader = (fnptr_glCreateShader)load_gl_fn(opengl32_dll_module, "glCreateShader");
        opengl_state->glShaderSource = (fnptr_glShaderSource)load_gl_fn(opengl32_dll_module, "glShaderSource");
        opengl_state->glCompileShader = (fnptr_glCompileShader)load_gl_fn(opengl32_dll_module, "glCompileShader");
        opengl_state->glGetShaderiv = (fnptr_glGetShaderiv)load_gl_fn(opengl32_dll_module, "glGetShaderiv");
        opengl_state->glGetShaderInfoLog = (fnptr_glGetShaderInfoLog)load_gl_fn(opengl32_dll_module, "glGetShaderInfoLog");
        opengl_state->glCreateProgram = (fnptr_glCreateProgram)load_gl_fn(opengl32_dll_module, "glCreateProgram");
        opengl_state->glAttachShader = (fnptr_glAttachShader)load_gl_fn(opengl32_dll_module, "glAttachShader");
        opengl_state->glLinkProgram = (fnptr_glLinkProgram)load_gl_fn(opengl32_dll_module, "glLinkProgram");
        opengl_state->glGetProgramiv = (fnptr_glGetProgramiv)load_gl_fn(opengl32_dll_module, "glGetProgramiv");
        opengl_state->glGetProgramInfoLog = (fnptr_glGetProgramInfoLog)load_gl_fn(opengl32_dll_module, "glGetProgramInfoLog");
        opengl_state->glDeleteShader = (fnptr_glDeleteShader)load_gl_fn(opengl32_dll_module, "glDeleteShader");
        opengl_state->glUseProgram = (fnptr_glUseProgram)load_gl_fn(opengl32_dll_module, "glUseProgram");
        opengl_state->glBufferSubData = (fnptr_glBufferSubData)load_gl_fn(opengl32_dll_module, "glBufferSubData");
        opengl_state->glDrawArrays = (fnptr_glDrawArrays)load_gl_fn(opengl32_dll_module, "glDrawArrays");
        opengl_state->glDrawArraysInstanced = (fnptr_glDrawArraysInstanced)load_gl_fn(opengl32_dll_module, "glDrawArraysInstanced");
        opengl_state->glDrawElementsInstanced = (fnptr_glDrawElementsInstanced)load_gl_fn(opengl32_dll_module, "glDrawElementsInstanced");
        opengl_state->glGetUniformLocation = (fnptr_glGetUniformLocation)load_gl_fn(opengl32_dll_module, "glGetUniformLocation");
        opengl_state->glUniformMatrix4fv = (fnptr_glUniformMatrix4fv)load_gl_fn(opengl32_dll_module, "glUniformMatrix4fv");
        opengl_state->glVertexAttribDivisor = (fnptr_glVertexAttribDivisor)load_gl_fn(opengl32_dll_module, "glVertexAttribDivisor");

        opengl_state->glGenTextures = (fnptr_glGenTextures)load_gl_fn(opengl32_dll_module, "glGenTextures");
        opengl_state->glBindTexture = (fnptr_glBindTexture)load_gl_fn(opengl32_dll_module, "glBindTexture");
        opengl_state->glTexParameteri = (fnptr_glTexParameteri)load_gl_fn(opengl32_dll_module, "glTexParameteri");
        opengl_state->glTexImage2D = (fnptr_glTexImage2D)load_gl_fn(opengl32_dll_module, "glTexImage2D");
        opengl_state->glTexSubImage2D = (fnptr_glTexSubImage2D)load_gl_fn(opengl32_dll_module, "glTexSubImage2D");

        glViewport(0, 0, (GLsizei)opengl_state->screen_width, (GLsizei)opengl_state->screen_height);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CCW);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);


        {
            // Vertex array object
            CALL_GL(glGenVertexArrays, 1, &opengl_state->vertex_array_object);
            CALL_GL(glBindVertexArray, opengl_state->vertex_array_object);

            // Index buffer object
            CALL_GL(glGenBuffers, 1, &opengl_state->index_buffer_object);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, opengl_state->index_buffer_object);
            CALL_GL(glBufferData, GL_ELEMENT_ARRAY_BUFFER, INDEX_ARRAY_BYTES, 0, GL_DYNAMIC_DRAW);

            u32_m attr_idx = 0;

            // Vertex buffer: vx
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_vx);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vx);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: vy
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_vy);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vy);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: vz
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_vz);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_vz);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: nx
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_nx);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nx);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: ny
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_ny);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_ny);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: nz
            CALL_GL(glGenBuffers, 1, &opengl_state->vertex_buffer_object_nz);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->vertex_buffer_object_nz);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Instanced vertex buffer: offset X
            CALL_GL(glGenBuffers, 1, &opengl_state->instanced_vertex_buffer_object_offset_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Instanced vertex buffer: offset Y
            CALL_GL(glGenBuffers, 1, &opengl_state->instanced_vertex_buffer_object_offset_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Instanced vertex buffer: offset Z
            CALL_GL(glGenBuffers, 1, &opengl_state->instanced_vertex_buffer_object_offset_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->instanced_vertex_buffer_object_offset_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

            // Shaders
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in float a_vx;\n"
                "layout (location = 1) in float a_vy;\n"
                "layout (location = 2) in float a_vz;\n"
                "layout (location = 3) in float a_nx;\n"
                "layout (location = 4) in float a_ny;\n"
                "layout (location = 5) in float a_nz;\n"
                "layout (location = 6) in float a_offset_x;\n"
                "layout (location = 7) in float a_offset_y;\n"
                "layout (location = 8) in float a_offset_z;\n"
                "uniform mat4 m_mvp;\n"
                "out vec3 v_normal;\n"
                "void main()\n"
                "{\n"
                "    vec3 pos = vec3(a_vx, a_vy, a_vz) + vec3(a_offset_x, a_offset_y, a_offset_z);\n"
                "    v_normal = vec3(a_nx, a_ny, a_nz);"
                "    gl_Position = m_mvp * vec4(pos, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "in vec3 v_normal;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    vec3 base_color = vec3(0.1f, 0.7f, 0.0f);\n"
                "    vec3 light_dir = -vec3(0.5f, 0.5f, 0.5f);\n"
                "    float inten = clamp(dot(light_dir, v_normal), 0.2f, 1.0f);\n"
                "    result_frag_color = vec4(base_color * inten, 1.0f);\n"
                //"    result_frag_color = vec4(0.05f, 0.5f, 0.0f, 1.0f);\n"
                //"    result_frag_color = vec4(0.05f, 0.5f, 0.0f, 1.0f);\n"
                "}\n"
                "\n";
            opengl_state->shader_program = make_shader_program(vertex_source, fragment_source);
        }


        {
            // Debug lines
            // Vertex array object
            CALL_GL(glGenVertexArrays, 1, &opengl_state->debug_line_vertex_array_object);
            CALL_GL(glBindVertexArray, opengl_state->debug_line_vertex_array_object);

            u32_m attr_idx = 0;

            // Vertex buffer: vertices
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_vertex_buffer_object_vertices);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_vertex_buffer_object_vertices);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, 6 * sizeof(f32), NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 3, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: start x
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: start y
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: start z
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_start_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end x
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end y
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: end z
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_end_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color r
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color g
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color b
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_line_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

            // Shaders
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in vec3 a_vertex;\n"
                "layout (location = 1) in float a_start_x;\n"
                "layout (location = 2) in float a_start_y;\n"
                "layout (location = 3) in float a_start_z;\n"
                "layout (location = 4) in float a_end_x;\n"
                "layout (location = 5) in float a_end_y;\n"
                "layout (location = 6) in float a_end_z;\n"
                "layout (location = 7) in float a_color_r;\n"
                "layout (location = 8) in float a_color_g;\n"
                "layout (location = 9) in float a_color_b;\n"
                "uniform mat4 m_mvp;\n"
                "out vec3 v_color;\n"
                "void main()\n"
                "{\n"
                "    vec3 start = vec3(a_start_x, a_start_y, a_start_z);\n"
                "    vec3 end = vec3(a_end_x, a_end_y, a_end_z);\n"
                "    vec3 pos = mix(start, end, a_vertex);\n"
                "    v_color = vec3(a_color_r, a_color_g, a_color_b);\n"
                "    gl_Position = m_mvp * vec4(pos, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "in vec3 v_color;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    result_frag_color = vec4(v_color, 1.0f);\n"
                "}\n"
                "\n";
            opengl_state->debug_line_shader_program = make_shader_program(vertex_source, fragment_source);
        }

        {
            // Debug spheres
            // Vertex array object
            CALL_GL(glGenVertexArrays, 1, &opengl_state->debug_sphere_vertex_array_object);
            CALL_GL(glBindVertexArray, opengl_state->debug_sphere_vertex_array_object);

            u32_m attr_idx = 0;

            // Vertex buffer: vertices
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_vertex_buffer_object_vertices);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_vertex_buffer_object_vertices);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, 12 * sizeof(f32), NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 3, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            // Vertex buffer: pos x
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_x);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_x);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: pos y
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_y);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_y);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: pos z
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_z);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_pos_z);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: radius
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_radius);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_radius);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color r
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_r);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color g
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_g);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            // Vertex buffer: color b
            CALL_GL(glGenBuffers, 1, &opengl_state->debug_sphere_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->debug_sphere_instanced_vertex_buffer_object_color_b);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, VERTEX_ARRAY_BYTES, NULL, GL_DYNAMIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 1, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            CALL_GL(glVertexAttribDivisor, attr_idx, 1);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);

            // Shaders
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in vec3 a_vertex;\n"
                "layout (location = 1) in float a_pos_x;\n"
                "layout (location = 2) in float a_pos_y;\n"
                "layout (location = 3) in float a_pos_z;\n"
                "layout (location = 4) in float a_radius;\n"
                "layout (location = 5) in float a_color_r;\n"
                "layout (location = 6) in float a_color_g;\n"
                "layout (location = 7) in float a_color_b;\n"
                "uniform mat4 m_cam;\n"
                "uniform mat4 m_proj;\n"
                "out vec3 v_color;\n"
                "out vec3 v_circle_center;\n"
                "out float v_circle_radius;\n"
                "out vec3 v_pos_worldspace;\n"
                "void main()\n"
                "{\n"
                "    v_color = vec3(a_color_r, a_color_g, a_color_b);\n"
                "    vec3 circle_center = vec3(a_pos_x, a_pos_y, a_pos_z);\n"
                "    vec3 vertex_pos_worldspace = \n"
                "        vec3(m_cam[0][0], m_cam[1][0], m_cam[2][0]) * a_vertex.x * a_radius + \n"
                "        vec3(m_cam[0][1], m_cam[1][1], m_cam[2][1]) * a_vertex.y * a_radius + \n"
                "        circle_center;    \n"
                "    v_circle_center = circle_center;\n"
                "    v_circle_radius = a_radius;\n"
                "    v_pos_worldspace = vertex_pos_worldspace;\n"
                "    gl_Position = m_proj * m_cam * vec4(vertex_pos_worldspace, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "in vec3 v_color;\n"
                "in vec3 v_circle_center;\n"
                "in float v_circle_radius;\n"
                "in vec3 v_pos_worldspace;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    vec3 d = v_pos_worldspace - v_circle_center;\n"
                "    if(dot(d, d) < v_circle_radius*v_circle_radius)\n"
                "    {\n"
                "        result_frag_color = vec4(v_color, 1.0f);\n"
                "    }\n"
                "    else\n"
                "    {\n"
                "        discard;\n"
                "    }\n"
                "    \n"
                "}\n"
                "\n";
            opengl_state->debug_sphere_shader_program = make_shader_program(vertex_source, fragment_source);
        }


        {
            CALL_GL(glGenTextures, 1, &opengl_state->fullscreen_quad_texture);
            CALL_GL(glBindTexture, GL_TEXTURE_2D, opengl_state->fullscreen_quad_texture);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            CALL_GL(glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            CALL_GL(glTexImage2D, GL_TEXTURE_2D,
                         0,
                         GL_RGBA,
                         opengl_state->screen_width,
                         opengl_state->screen_height,
                         0,
                         GL_RGBA,
                         GL_UNSIGNED_BYTE,
                         NULL);
            const char* vertex_source =
                "#version 330 core\n"
                "layout (location = 0) in vec2 a_pos;\n"
                "out vec2 v_uv;\n"
                "void main()\n"
                "{\n"
                "    v_uv = a_pos * 0.5f + vec2(0.5f, 0.5f);\n"
                "    gl_Position = vec4(a_pos, 0.0f, 1.0f);\n"
                "}\n"
                "\n";
            const char* fragment_source =
                "#version 330 core\n"
                "uniform sampler2D tex_sampler;\n"
                "in vec2 v_uv;\n"
                "out vec4 result_frag_color;\n"
                "void main()\n"
                "{\n"
                "    result_frag_color = texture(tex_sampler, v_uv);\n"
                "}\n"
                "\n";
            opengl_state->textured_quad_shader_program = make_shader_program(vertex_source, fragment_source);

            CALL_GL(glGenVertexArrays, 1, &opengl_state->fullscreen_quad_vertex_array_object);
            CALL_GL(glBindVertexArray, opengl_state->fullscreen_quad_vertex_array_object);

            u32_m attr_idx = 0;

            // Vertex buffer: vertices
            CALL_GL(glGenBuffers, 1, &opengl_state->fullscreen_quad_vertex_buffer_object_pos);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, opengl_state->fullscreen_quad_vertex_buffer_object_pos);
            CALL_GL(glBufferData, GL_ARRAY_BUFFER, 6 * 2 * sizeof(f32), NULL, GL_STATIC_DRAW);
            CALL_GL(glVertexAttribPointer, attr_idx, 2, GL_FLOAT, GL_FALSE, 0, 0);
            CALL_GL(glEnableVertexAttribArray, attr_idx);
            attr_idx++;

            CALL_GL(glBindVertexArray, 0);
            CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
        }
    }

    // Init engine input state.
    // Dependent on window being created.
    {
        struct InputState* input_state = &g_platform->input_state;
        POINT p;
        const BOOL get_cursor_pos_succes = GetCursorPos(&p);
        ASSERT(get_cursor_pos_succes, "Failed GetCursorPos");
        const BOOL screen_to_client_success = ScreenToClient(platform->hwnd, &p);
        ASSERT(screen_to_client_success, "Failed ScreenToClient");
        input_state->mouse_screen_pos_x = (s32_m)p.x;
        input_state->mouse_screen_pos_y = (s32_m)p.y;
        input_state->mouse_screen_dx = 0;
        input_state->mouse_screen_dy = 0;
        input_state->fps_mode = 0;
    }

    // Init time.
    {
        LARGE_INTEGER clock_freq;
        QueryPerformanceFrequency(&clock_freq);
        platform->clock_freq = clock_freq.QuadPart;
    }


    u64 frame_memory_arena_bytes = MB(500);
    struct MemoryArena frame_memory_arena = memory_arena_init(MEMORY_ARENA_ALLOCATE(&main_memory_arena, frame_memory_arena_bytes), frame_memory_arena_bytes);
    struct GameState* game_state_buffer = (struct GameState*)MEMORY_ARENA_ALLOCATE(&main_memory_arena, 2 * sizeof(struct GameState));
    struct GameState* game_state_A = game_state_buffer;
    struct GameState* game_state_B = game_state_buffer + 1;
    reset_game_state(game_state_A);
    reset_game_state(game_state_B);

    // Main loop
    s64 engine_start_timestamp_us = get_timestamp_us();
    s64_m last_frame_timestamp_us = engine_start_timestamp_us;
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


        // Main loop time control.
#define ENGINE_FRAME_DURATION_US 8333LL
        {
            s64 time_since_last_frame_us = get_timestamp_us() - last_frame_timestamp_us;
            if(time_since_last_frame_us >= ENGINE_FRAME_DURATION_US)
            {
                last_frame_timestamp_us += ENGINE_FRAME_DURATION_US;
            }
            else
            {
                continue;
            }
        }

        // Read input
        {
            struct InputState* input_state = &platform->input_state;
            s32_m last_mouse_screen_pos_x = input_state->mouse_screen_pos_x;
            s32_m last_mouse_screen_pos_y = input_state->mouse_screen_pos_y;
            POINT p;
            const BOOL get_cursor_pos_succes = GetCursorPos(&p);
            ASSERT(get_cursor_pos_succes, "Failed GetCursorPos");
            input_state->mouse_screen_pos_x = (s32_m)p.x;
            input_state->mouse_screen_pos_y = (s32_m)p.y;

            if(is_key_toggled_down(KB_E))
            {
                input_state->fps_mode = !input_state->fps_mode;
            }

            if(input_state->fps_mode)
            {
                show_cursor(0);
                RECT clip_rect;
                const BOOL get_client_rect_success = GetClientRect(platform->hwnd, &clip_rect);
                ASSERT(get_client_rect_success, "GetClientRect failed.");
                s32 tx = (clip_rect.left + clip_rect.right) / 2;
                s32 ty = (clip_rect.bottom + clip_rect.top) / 2;

                POINT client_to_screen_point = {
                    .x = tx,
                    .y = ty,
                };
                const BOOL client_to_screen_success = ClientToScreen(platform->hwnd, &client_to_screen_point);
                ASSERT(client_to_screen_success, "ClientToScreen failed.");

                input_state->mouse_screen_dx = input_state->mouse_screen_pos_x - last_mouse_screen_pos_x;
                input_state->mouse_screen_dy = input_state->mouse_screen_pos_y - last_mouse_screen_pos_y;

                const BOOL set_cursor_pos_success = SetCursorPos(client_to_screen_point.x, client_to_screen_point.y);
                ASSERT(set_cursor_pos_success, "SetCursorPos failed.");

                input_state->mouse_screen_pos_x = client_to_screen_point.x;
                input_state->mouse_screen_pos_y = client_to_screen_point.y;
            }
            else
            {
                show_cursor(1);
                input_state->mouse_screen_dx = input_state->mouse_screen_pos_x - last_mouse_screen_pos_x;
                input_state->mouse_screen_dy = input_state->mouse_screen_pos_y - last_mouse_screen_pos_y;
            }
        }

        reset_game_state(game_state_A);

        do_one_frame(
            game_state_A,
            game_state_B,
            &frame_memory_arena
        );

        draw_game_state(game_state_A, &frame_memory_arena);

        // Rotate game state double buffer.
        {
            struct GameState* tmp = game_state_A;
            game_state_A = game_state_B;
            game_state_B = tmp;
        }

        memory_arena_reset(&frame_memory_arena);

        {
            struct InputState* input_state = &g_platform->input_state;
            memcpy(input_state->last_key, input_state->key, sizeof(input_state->last_key));
            memcpy(input_state->last_mouse_key, input_state->mouse_key, sizeof(input_state->last_mouse_key));
        }

        const HDC dc = GetDC(platform->hwnd);
        BOOL swap_buffers_success = SwapBuffers(dc);
        ASSERT(swap_buffers_success, "SwapBuffers failed.");

        glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    CloseHandle(g_platform->log_file);

    ExitProcess(0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

