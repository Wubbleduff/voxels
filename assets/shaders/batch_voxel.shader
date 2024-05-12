
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
    uint b_color[10*1024];
};

uniform mat4 m_view;
uniform mat4 m_proj;

out flat vec4 model_position;
out flat vec4 world_position;
out flat vec4 view_position;
out flat vec4 model_color;
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

    uint color_tmp = b_color[gl_InstanceID];
    uint r = (color_tmp & 0xFF000000) >> 24;
    uint g = (color_tmp & 0x00FF0000) >> 16;
    uint b = (color_tmp & 0x0000FF00) >> 8;
    model_color = vec4(r / 256.0f, g / 256.0f, b / 256.0f, 1.0f);

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

    float inten = dot(normalize(vec4(1.0f, 1.0f, 1.0f, 0.0f)), n);

    frag_color = vec4(model_color.x, model_color.y, model_color.z, 1.0f);
};

