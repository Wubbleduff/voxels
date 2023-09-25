
@vertex



#version 440 core
layout (location = 0) in float a_x;
layout (location = 1) in float a_y;
layout (location = 2) in float a_z;
layout (location = 3) in float n_x;
layout (location = 4) in float n_y;
layout (location = 5) in float n_z;

layout(std430, binding = 0) buffer voxel_data
{
    float b_px[10*1024];
    float b_py[10*1024];
    float b_pz[10*1024];
    float b_scale[10*1024];
    int b_color[10*1024];
};

uniform mat4 m_view;
uniform mat4 m_proj;

out vec4 model_position;
out vec4 world_position;
out vec4 view_position;
out vec4 model_color;
out vec4 normal;

void main()
{
    float offset_x = b_px[gl_InstanceID];
    float offset_y = b_py[gl_InstanceID];
    float offset_z = b_pz[gl_InstanceID];
    float scale = b_scale[gl_InstanceID];

    model_position = vec4(a_x, a_y, a_z, 1.0f);
    world_position = vec4(a_x*scale + offset_x, a_y*scale + offset_y, a_z*scale + offset_z, 1.0f);
    view_position = m_view * world_position;

    int color_tmp = b_color[gl_InstanceID];
    int r = (color_tmp & 0xFF000000) >> 24;
    int g = (color_tmp & 0x00FF0000) >> 16;
    int b = (color_tmp & 0x0000FF00) >> 8;
    int a = (color_tmp & 0x000000FF) >> 0;
    model_color = vec4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);

    normal = vec4(n_x, n_y, n_z, 1.0f);
    gl_Position = m_proj * view_position;
};



@fragment



#version 440 core


uniform float u_view_dist;

in vec4 model_position;
in vec4 world_position;
in vec4 view_position;
in vec4 model_color;
in vec4 normal;

out vec4 frag_color;

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
    vec4 n;
    if(abs(normal.x) > abs(normal.y) && abs(normal.x) > abs(normal.z))
    {
        n = vec4(1.0f, 0.0f, 0.0f, 0.0f);
        n *= sign(normal.x);
    }
    else if(abs(normal.y) > abs(normal.x) && abs(normal.y) > abs(normal.z))
    {
        n = vec4(0.0f, 1.0f, 0.0f, 0.0f);
        n *= sign(normal.y);
    }
    else if(abs(normal.z) > abs(normal.x) && abs(normal.z) > abs(normal.y))
    {
        n = vec4(0.0f, 0.0f, 1.0f, 0.0f);
        n *= sign(normal.z);
    }

    //vec4 n = normalize(normal);
    float inten = dot(normalize(vec4(1.0f, 1.0f, 1.0f, 0.0f)), n);
    //inten = max(inten, 0.2f);

    frag_color = model_color * inten;
    //frag_color = model_color;
};

