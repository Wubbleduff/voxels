
@vertex



#version 440 core
layout (location = 0) in float a_x;
layout (location = 1) in float a_y;
layout (location = 2) in float a_z;

layout(std430, binding = 0) buffer voxel_data
{
    float px[65536];
    float py[65536];
    float pz[65536];
};

uniform mat4 vp;

out vec4 model_position;
out vec4 world_position;

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
    gl_Position = vp * world_position;
};



@fragment



#version 440 core

in vec4 model_position;
in vec4 world_position;

out vec4 frag_color;

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
    vec4 light_pos = vec4(0.0f, 5.0f, 0.0f, 1.0f);
    float inten = dot(normalize(light_pos - world_position), normal);
    frag_color = vec4(inten, inten, inten, 1.0f);
};

