layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 texcoord;

uniform mat4 view_projection;

out vec3 PS_IN_WorldPos;

void main(void)
{
    PS_IN_WorldPos = position;
    gl_Position = view_projection * vec4(position, 1.0f);
}
