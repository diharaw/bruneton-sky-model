layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 texcoords;

out vec3 PS_IN_TexCoord;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    PS_IN_TexCoord = position;
    vec4 pos = projection * view * vec4(position, 1.0);
    gl_Position = pos.xyww;
}  