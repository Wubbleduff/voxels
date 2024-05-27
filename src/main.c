
#include "common.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <immintrin.h>

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
// Math

#define H_PI 1.57079637f
#define PI 3.14159274f
#define TAU 6.28318548f
#define INV_PI 0.318309873f
#define INV_H_PI 0.636620f
#define INV_TAU 0.159154937f

typedef struct
{
    union
    {
        _Alignas(32) f32_m m[16];
        __m128 v[4];
    };
} mtx4x4_m;
typedef const mtx4x4_m mtx4x4;

typedef struct
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
    };
} v2_m;
typedef const v2_m v2;

typedef struct
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
    };
} v3_m;
typedef const v3_m v3;

typedef struct
{
    union
    {
        _Alignas(16) f32_m m[4];
        __m128 v;
    };
} v4_m;
typedef const v4_m v4;

MAYBE_UNUSED static inline f32 v3_dot(v3 a, v3 b)
{
    // | a.x | a.y | a.z | - |
    // | b.x | b.y | b.z | - |

    v3_m r;

    //     00        01        10
    // | a.x*b.x | a.y*b.y | a.z*b.z | - |
    __m128 p = _mm_mul_ps(a.v, b.v);

    // | a.x*b.x | a.y*b.y | a.z*b.z | - |
    // +
    // | a.y*b.y | a.z*b.z | a.x*b.x | - |
    r.v = _mm_add_ps(p, _mm_shuffle_ps(p, p, 0b11001001));

    // | a.x*b.x + a.y*b.y | a.y*b.y + a.z*b.z | a.z*b.z + a.x*b.x | - |
    // +
    // |           a.z*b.z |           a.x*b.x |           a.y*b.y | - |
    r.v = _mm_add_ps(r.v, _mm_shuffle_ps(p, p, 0b11010010));

    return r.m[0];
}

MAYBE_UNUSED static inline v3 v3_cross(v3 a, v3 b)
{
    v3_m r;
    r.v = _mm_fmsub_ps(
        _mm_shuffle_ps(a.v, a.v, 0b11001001),
        _mm_shuffle_ps(b.v, b.v, 0b11010010),
        _mm_mul_ps(
            _mm_shuffle_ps(a.v, a.v, 0b11010010),
            _mm_shuffle_ps(b.v, b.v, 0b11001001)
        )
    );
    return r;
}

MAYBE_UNUSED static inline v3 v3_add(v3 a, v3 b)
{
    v3_m r;
    r.v = _mm_add_ps(a.v, b.v);
    return r;
}

MAYBE_UNUSED static inline v3 v3_scale(v3 a, f32 b)
{
    v3_m r;
    r.v = _mm_mul_ps(a.v, _mm_set1_ps(b));
    return r;
}

MAYBE_UNUSED static inline v3 v3_normalize(v3 a)
{
    // TODO(mfritz) No need to go to scalar here.
    v3_m r = a;
    f32 len_sqr = v3_dot(a, a);
    f32 len = _mm_cvtss_f32(_mm_sqrt_ss(_mm_set1_ps(len_sqr)));
    r = v3_scale(r, 1.0f / len);
    return r;
}

MAYBE_UNUSED static inline __m256 approx_sin8(__m256 x)
{
    const __m256 pi       = _mm256_set1_ps(PI);
    const __m256 h_pi     = _mm256_set1_ps(H_PI);
    const __m256 inv_pi   = _mm256_set1_ps(INV_PI);
    const __m256 inv_h_pi = _mm256_set1_ps(INV_H_PI);
    const __m256 inv_tau = _mm256_set1_ps(INV_TAU);
    
    // Range reduce to [0, 2*pi] to check if we need to negate result.
    __m256 sign;
    {
        __m256 xx = _mm256_mul_ps(x, inv_tau);
        xx = _mm256_sub_ps(xx, _mm256_round_ps(xx, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
        __m256 mask = _mm256_cmp_ps(xx, _mm256_set1_ps(0.5f), _CMP_NLT_UQ);
        sign = _mm256_and_ps(mask, _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));
    }
    
    // From now on x is positive
    x = _mm256_and_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)));
    
    // Range reduce to [0, pi] to check if we need to flip X.
    __m256 flip_mask;
    {
        __m256 xx = _mm256_mul_ps(x, inv_pi);
        xx = _mm256_sub_ps(xx, _mm256_round_ps(xx, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
        flip_mask = _mm256_cmp_ps(xx, _mm256_set1_ps(0.5f), _CMP_NLT_UQ);
    }
    
    // Range reduce to [0, pi/2] for evaluation
    x = _mm256_mul_ps(x, inv_h_pi);
    x = _mm256_sub_ps(x, _mm256_round_ps(x, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    x = _mm256_mul_ps(x,  h_pi);
    
    __m256 flip_x = _mm256_fmadd_ps(x, _mm256_set1_ps(-1.0f), h_pi);
    x = _mm256_blendv_ps(x, flip_x, flip_mask);
    
    // x = x - x^3 / 6
    __m256 eval0 =  _mm256_fnmadd_ps(_mm256_mul_ps(_mm256_mul_ps(x, x), x), _mm256_set1_ps(1.0f / 6.0f), x);
    
    // x = 1/384*(pi - 2*x)^4 - 1/8*(pi - 2x)^2 + 1
    __m256 t = _mm256_fnmadd_ps(x, _mm256_set1_ps(2.0f), pi);
    __m256 t2 = _mm256_mul_ps(t, t);
    __m256 t4 = _mm256_mul_ps(t2, t2);
    __m256 eval1 = _mm256_mul_ps(_mm256_set1_ps(1.0f/384.0f), t4);
    eval1 = _mm256_fnmadd_ps(_mm256_set1_ps(1.0f/8.0f), t2, eval1);
    eval1 = _mm256_add_ps(eval1, _mm256_set1_ps(1.0f));
    
    __m256 eval_blend_mask = _mm256_cmp_ps(x, _mm256_set1_ps(0.6403f), _CMP_LT_OQ);
    __m256 eval = _mm256_blendv_ps(eval1, eval0, eval_blend_mask);
    
    eval = _mm256_or_ps(eval, sign);
    
    return eval;
}

// 4x4 matrix multiply : r = a * b.
// NOTE: Matrices are assumed to be row-major.
MAYBE_UNUSED static inline void mtx4x4_mul(mtx4x4_m* r, mtx4x4* a, mtx4x4* b)
{
    /*
     *
     * | a00  a01  a02  a03 |   | b00  b01  b02  b03 |   
     * | a10  a11  a12  a13 | * | b10  b11  b12  b13 | = 
     * | a20  a21  a22  a23 |   | b20  b21  b22  b23 |   
     * | a30  a31  a32  a33 |   | b30  b31  b32  b33 |   
     *
     * | a00*b00 + a01*b10 + a02*b20 + a03*b30    a00*b01 + a01*b11 + a02*b21 + a03*b31    a00*b02 + a01*b12 + a02*b22 + a03*b32    a00*b03 + a01*b13 + a02*b23 + a03*b33 |
     * | a10*b00 + a11*b10 + a12*b20 + a13*b30    a10*b01 + a11*b11 + a12*b21 + a13*b31    a10*b02 + a11*b12 + a12*b22 + a13*b32    a10*b03 + a11*b13 + a12*b23 + a13*b33 |
     * | a20*b00 + a21*b10 + a22*b20 + a23*b30    a20*b01 + a21*b11 + a22*b21 + a23*b31    a20*b02 + a21*b12 + a22*b22 + a23*b32    a20*b03 + a21*b13 + a22*b23 + a23*b33 |
     * | a30*b00 + a31*b10 + a32*b20 + a33*b30    a30*b01 + a31*b11 + a32*b21 + a33*b31    a30*b02 + a31*b12 + a32*b22 + a33*b32    a30*b03 + a31*b13 + a32*b23 + a33*b33 |
     *
     */

    /*
    Reference
    for(u64_m i = 0; i < 4; i++)
    {
        for(u64_m j = 0; j < 4; j++)
        {
            r[i*4 + j] =
                a[i*4 + 0] * b[0*4 + j] + 
                a[i*4 + 1] * b[1*4 + j] + 
                a[i*4 + 2] * b[2*4 + j] + 
                a[i*4 + 3] * b[3*4 + j];
        }
    }
    */

    for(u64_m row = 0; row < 4; row++)
    {
        const __m128 v = _mm_fmadd_ps(
            _mm_broadcast_ss(a->m + row*4 + 0),
            _mm_loadu_ps(b->m + 0*4 + 0),
            _mm_fmadd_ps(
                _mm_broadcast_ss(a->m + row*4 + 1),
                _mm_loadu_ps(b->m + 1*4 + 0),
                _mm_fmadd_ps(
                    _mm_broadcast_ss(a->m + row*4 + 2),
                    _mm_loadu_ps(b->m + 2*4 + 0),
                    _mm_mul_ps(
                        _mm_broadcast_ss(a->m + row*4 + 3),
                        _mm_loadu_ps(b->m + 3*4 + 0)
                    )
                )
            )
        );
        _mm_storeu_ps(r->m + row*4, v);
    }
}

#if 1
MAYBE_UNUSED static inline void make_x_axis_rotation_mtx(mtx4x4_m* r, f32 turns)
{
    __m256 a_v = approx_sin8(_mm256_setr_ps(TAU * turns, TAU * turns + H_PI, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
    _Alignas(32) f32_m a[8];
    _mm256_store_ps(a, a_v);
    f32 sin_a = a[0];
    f32 cos_a = a[1];

    r->m[0] = 1.0f;   r->m[1] = 0.0f;   r->m[2] = 0.0f;    r->m[3] = 0.0f;
    r->m[4] = 0.0f;   r->m[5] = cos_a;  r->m[6] = -sin_a;  r->m[7] = 0.0f;
    r->m[8] = 0.0f;   r->m[9] = sin_a;  r->m[10] = cos_a;  r->m[11] = 0.0f;
    r->m[12] = 0.0f;  r->m[13] = 0.0f;  r->m[14] = 0.0f;   r->m[15] = 1.0f;
}

MAYBE_UNUSED static inline void make_y_axis_rotation_mtx(mtx4x4_m* r, float turns)
{
    __m256 a_v = approx_sin8(_mm256_setr_ps(TAU * turns, TAU * turns + H_PI, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f));
    _Alignas(32) f32_m a[8];
    _mm256_store_ps(a, a_v);
    f32 sin_a = a[0];
    f32 cos_a = a[1];

    r->m[0] = cos_a;   r->m[1] = 0.0f;   r->m[2] = sin_a;   r->m[3] = 0.0f;
    r->m[4] = 0.0f;    r->m[5] = 1.0f;   r->m[6] = 0.0f;    r->m[7] = 0.0f;
    r->m[8] = -sin_a;  r->m[9] = 0.0f;   r->m[10] = cos_a;  r->m[11] = 0.0f;
    r->m[12] = 0.0f;   r->m[13] = 0.0f;  r->m[14] = 0.0f;   r->m[15] = 1.0f;
}

MAYBE_UNUSED static inline void make_translation_mtx(mtx4x4_m* r, v3 v)
{
    r->m[0]  = 1.0f;   r->m[1] = 0.0f;    r->m[2] = 0.0f;   r->m[3] = v.m[0];
    r->m[4]  = 0.0f;   r->m[5] = 1.0f;    r->m[6] = 0.0f;   r->m[7] = v.m[1];
    r->m[8]  = 0.0f;   r->m[9] = 0.0f;   r->m[10] = 1.0f;  r->m[11] = v.m[2];
    r->m[12] = 0.0f;  r->m[13] = 0.0f;   r->m[14] = 0.0f;  r->m[15] = 1.0f;
}

#endif

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
typedef GLint (*fnptr_glGetUniformLocation)(GLuint program, const GLchar *name);
typedef void (*fnptr_glUniformMatrix4fv)(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);


#define VERTEX_ARRAY_BYTES MB(1)
#define INDEX_ARRAY_BYTES MB(1)

struct OpenGLState
{
    HWND hwnd;
    u64_m screen_width;
    u64_m screen_height;

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
    fnptr_glGetUniformLocation glGetUniformLocation;
    fnptr_glUniformMatrix4fv glUniformMatrix4fv;
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
    type MACRO_r = g_opengl_state->fn(__VA_ARGS__); \
    const GLuint err = g_opengl_state->glGetError(); \
    g_opengl_state->last_gl_error = err; \
    ASSERT(err == 0, #fn " failed"); \
    *(out) = MACRO_r; \
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
    // https://learn.microsoft.com/en-us/windows/win32/inputdev/virtual-key-codes
    KB_ESCAPE = VK_ESCAPE,
    KB_SPACE = VK_SPACE,
    KB_LCTRL = VK_CONTROL,
    KB_0 = 0x30,    
    KB_1 = 0x31,    
    KB_2 = 0x32,    
    KB_3 = 0x33,    
    KB_4 = 0x34,    
    KB_5 = 0x35,    
    KB_6 = 0x36,    
    KB_7 = 0x37,    
    KB_8 = 0x38,    
    KB_9 = 0x39,    
    KB_A = 0x41,    
    KB_B = 0x42,    
    KB_C = 0x43,    
    KB_D = 0x44,    
    KB_E = 0x45,    
    KB_F = 0x46,    
    KB_G = 0x47,    
    KB_H = 0x48,    
    KB_I = 0x49,    
    KB_J = 0x4A,    
    KB_K = 0x4B,    
    KB_L = 0x4C,    
    KB_M = 0x4D,    
    KB_N = 0x4E,    
    KB_O = 0x4F,    
    KB_P = 0x50,    
    KB_Q = 0x51,    
    KB_R = 0x52,    
    KB_S = 0x53,    
    KB_T = 0x54,    
    KB_U = 0x55,    
    KB_V = 0x56,    
    KB_W = 0x57,    
    KB_X = 0x58,    
    KB_Y = 0x59,    
    KB_Z = 0x5A,    
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

        case WM_SIZE:
        {
            const UINT width = LOWORD(lParam);
            const UINT height = HIWORD(lParam);
            g_opengl_state->screen_width = width;
            g_opengl_state->screen_height = height;
            glViewport(0, 0, width, height);
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
    memset(g_input_state->key, 0, sizeof(g_input_state->key));
    memset(g_input_state->mouse_key, 0, sizeof(g_input_state->mouse_key));

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

        RECT client_rect;
        const BOOL get_client_rect_success = GetClientRect(g_opengl_state->hwnd, &client_rect);
        ASSERT(get_client_rect_success, "GetClientRect failed.");
        g_opengl_state->screen_width = client_rect.right;
        g_opengl_state->screen_height = client_rect.bottom;

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
        g_opengl_state->glGetUniformLocation = (fnptr_glGetUniformLocation)load_gl_fn(opengl32_dll_module, "glGetUniformLocation");
        g_opengl_state->glUniformMatrix4fv = (fnptr_glUniformMatrix4fv)load_gl_fn(opengl32_dll_module, "glUniformMatrix4fv");

        glViewport(0, 0, (GLsizei)g_opengl_state->screen_width, (GLsizei)g_opengl_state->screen_height);

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
            "uniform mat4 m_proj;\n"
            "void main()\n"
            "{\n"
            "    vec3 pos = vec3(a_x, a_y, a_z);\n"
            "    gl_Position = m_proj * vec4(pos, 1.0f);\n"
            "}\n"
            "\n";
        const char* fragment_source =
            "#version 440 core\n"
            "out vec4 result_frag_color;\n"
            "void main()\n"
            "{\n"
            "    result_frag_color = vec4(0.1f, 0.0f, 0.4f, 1.0f);\n"
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

        f32 vx_offset = 0.0f;
        f32 vy_offset = 0.0f;
        f32 vz_offset = -3.0f;
        f32 vx[] = { vx_offset + 0.0f, vx_offset + 1.0f, vx_offset + 0.0f, vx_offset + 1.0f, vx_offset + 0.0f, vx_offset + 1.0f, vx_offset + 0.0f, vx_offset + 1.0f };
        f32 vy[] = { vy_offset + 0.0f, vy_offset + 0.0f, vy_offset + 1.0f, vy_offset + 1.0f, vy_offset + 0.0f, vy_offset + 0.0f, vy_offset + 1.0f, vy_offset + 1.0f };
        f32 vz[] = { vz_offset + 0.0f, vz_offset + 0.0f, vz_offset + 0.0f, vz_offset + 0.0f, vz_offset + 1.0f, vz_offset + 1.0f, vz_offset + 1.0f, vz_offset + 1.0f };
        u32 indices[] = 
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

        glClearColor(0.0f, 161.0f/255.0f, 201.0f/255.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        CALL_GL(glUseProgram, g_opengl_state->shader_program);


        {
            // Camera controls
            static f32_m cam_pitch_turns = 0.0f;
            cam_pitch_turns += (float)is_key_down(KB_I) * 0.001f;
            cam_pitch_turns -= (float)is_key_down(KB_K) * 0.001f;

            static f32_m cam_yaw_turns = 0.0f;
            cam_yaw_turns += (float)is_key_down(KB_J) * 0.001f;
            cam_yaw_turns -= (float)is_key_down(KB_L) * 0.001f;

            static v3_m cam_pos = {0};

            mtx4x4_m y_rot_mtx;
            make_y_axis_rotation_mtx(&y_rot_mtx, -cam_yaw_turns);

            mtx4x4_m x_rot_mtx;
            make_x_axis_rotation_mtx(&x_rot_mtx, -cam_pitch_turns);

            mtx4x4_m translation_mtx;
            make_translation_mtx(&translation_mtx, v3_scale(cam_pos, -1.0f));

            mtx4x4_m world_to_cam_mtx_temp;
            mtx4x4_mul(&world_to_cam_mtx_temp, &y_rot_mtx, &translation_mtx);

            mtx4x4_m world_to_cam_mtx;
            mtx4x4_mul(&world_to_cam_mtx, &x_rot_mtx, &world_to_cam_mtx_temp);

            v3 I_cam = {.m={world_to_cam_mtx.m[0*4 + 0], world_to_cam_mtx.m[0*4 + 1], world_to_cam_mtx.m[0*4 + 2]}};
            v3 J_cam = {.m={world_to_cam_mtx.m[1*4 + 0], world_to_cam_mtx.m[1*4 + 1], world_to_cam_mtx.m[1*4 + 2]}};
            v3 K_cam = {.m={world_to_cam_mtx.m[2*4 + 0], world_to_cam_mtx.m[2*4 + 1], world_to_cam_mtx.m[2*4 + 2]}};

            f32 speed = 0.1f;
            cam_pos = v3_add(cam_pos, v3_scale(K_cam, -speed * (float)is_key_down(KB_W)));
            cam_pos = v3_add(cam_pos, v3_scale(K_cam,  speed * (float)is_key_down(KB_S)));

            cam_pos = v3_add(cam_pos, v3_scale(I_cam,  speed * (float)is_key_down(KB_D)));
            cam_pos = v3_add(cam_pos, v3_scale(I_cam, -speed * (float)is_key_down(KB_A)));
            
            cam_pos = v3_add(cam_pos, v3_scale(J_cam,  speed * (float)is_key_down(KB_SPACE)));
            cam_pos = v3_add(cam_pos, v3_scale(J_cam, -speed * (float)is_key_down(KB_LCTRL)));

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

            f32 aspect_ratio = (float)g_opengl_state->screen_height / (float)g_opengl_state->screen_width;

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

            mtx4x4_m cam_proj_mtx;
            mtx4x4_mul(&cam_proj_mtx, &proj_mtx, &world_to_cam_mtx);

            GLint loc;
            CALL_GL_RET(&loc, GLint, glGetUniformLocation, g_opengl_state->shader_program, "m_proj");
            CALL_GL(glUniformMatrix4fv, loc, 1, 1, &cam_proj_mtx.m[0]);
            ASSERT(loc != -1, "Failed to bind uniform.");
        }


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

        CALL_GL(glBindVertexArray, 0);
        CALL_GL(glBindBuffer, GL_ARRAY_BUFFER, 0);
        CALL_GL(glBindBuffer, GL_ELEMENT_ARRAY_BUFFER, 0);

        const HDC dc = GetDC(g_opengl_state->hwnd);
        BOOL swap_buffers_success = SwapBuffers(dc);
        ASSERT(swap_buffers_success, "SwapBuffers failed.");
    }

    CloseHandle(g_log_file);

    ExitProcess(0);
}

