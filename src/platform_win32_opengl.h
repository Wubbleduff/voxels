
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <GL/GL.h>
#include "wglext.h"

#pragma once

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
typedef void (*fnptr_glDrawElements)(GLenum mode, GLsizei count, GLenum type, const void * indices);
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
    u32 screen_width;
    u32 screen_height;

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
    fnptr_glDrawElements glDrawElements;
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

INTERNAL const f32 s_cube_mesh_vx[] = { -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f, -0.5f,  0.5f };
INTERNAL const f32 s_cube_mesh_vy[] = { -0.5f, -0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f };
INTERNAL const f32 s_cube_mesh_vz[] = { -0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f };
INTERNAL const u32 s_cube_mesh_indices[] = 
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
    u32 key[NUM_KEYBOARD_KEYS];
    u32 mouse_key[MAX_MOUSE_KEYS];

    u32 last_key[NUM_KEYBOARD_KEYS];
    u32 last_mouse_key[MAX_MOUSE_KEYS];

    u32 fps_mode; 
    s32 mouse_screen_pos_x;
    s32 mouse_screen_pos_y;
    s32 mouse_screen_dx;
    s32 mouse_screen_dy;
};

// Maps Win32 VK code to KeyboardKey enum.
INTERNAL const u32 win32_vk_code_to_keyboard_key[] = {
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


struct Platform_Win32
{
    HANDLE log_file;
    HWND hwnd;
    u64 clock_freq;

    struct OpenGLState opengl_state;
    struct InputState input_state;
};


