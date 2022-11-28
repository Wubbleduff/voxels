
@vertex



#version 440 core
layout (location = 0) in float a_x;
layout (location = 1) in float a_y;
layout (location = 2) in float a_z;

layout(std430, binding = 0) buffer voxel_data
{
    float px[1024*1024];
    float py[1024*1024];
    float pz[1024*1024];
    int color[1024*1024];
};

uniform mat4 m_view;
uniform mat4 m_proj;

out vec4 model_position;
out vec4 world_position;
out vec4 view_position;
out vec4 model_color;

void main()
{
    float offset_x = px[gl_InstanceID];
    float offset_y = py[gl_InstanceID];
    float offset_z = pz[gl_InstanceID];

    // TODO speed
    mat4 translate = mat4(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        offset_x, offset_y, offset_z, 1.0f
    );

    mat4 world_m_object = translate;

    model_position = vec4(a_x, a_y, a_z, 1.0f);
    world_position = world_m_object * vec4(a_x, a_y, a_z, 1.0f);
    view_position = m_view * world_position;






    int color_tmp = color[gl_InstanceID];

    int r = (color_tmp & 0xFF000000) >> 24;
    int g = (color_tmp & 0x00FF0000) >> 16;
    int b = (color_tmp & 0x0000FF00) >> 8;
    int a = (color_tmp & 0x000000FF) >> 0;
    model_color = vec4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);

    gl_Position = m_proj * view_position;
};



@fragment



#version 440 core


uniform float u_view_dist;

in vec4 model_position;
in vec4 world_position;
in vec4 view_position;
in vec4 model_color;

out vec4 frag_color;

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
    vec4 normal;
    if(abs(model_position.x) > abs(model_position.y) && abs(model_position.x) > abs(model_position.z))
    {
        normal = vec4(1.0f, 0.0f, 0.0f, 0.0f);
        normal *= sign(model_position.x);
    }
    else if(abs(model_position.y) > abs(model_position.x) && abs(model_position.y) > abs(model_position.z))
    {
        normal = vec4(0.0f, 1.0f, 0.0f, 0.0f);
        normal *= sign(model_position.y);
    }
    else if(abs(model_position.z) > abs(model_position.x) && abs(model_position.z) > abs(model_position.y))
    {
        normal = vec4(0.0f, 0.0f, 1.0f, 0.0f);
        normal *= sign(model_position.z);
    }
    float inten = dot(vec4(0.0f, 1.0f, 0.0f, 0.0f), normal);
    inten = max(inten, 0.2f);

    frag_color = model_color * inten;
};

