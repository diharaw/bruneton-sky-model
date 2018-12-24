layout (location = 0) in vec3 VS_IN_Position;
layout (location = 1) in vec2 VS_IN_Texcoord;
layout (location = 2) in vec3 VS_IN_Normal;
layout (location = 3) in vec3 VS_IN_Tangent;
layout (location = 4) in vec3 VS_IN_Bitangent;

layout (std140) uniform u_GlobalUBO
{ 
    mat4 view;
    mat4 projection;
	mat4 inv_view;
	mat4 inv_projection;
	vec4 camera_pos;
};

uniform mat4 model;

out vec2 PS_IN_TexCoord;
out vec3 PS_IN_FragPos;
out vec3 PS_IN_Normal;

void main(void)
{
    vec4 world_pos = model * vec4(VS_IN_Position, 1.0f);
    PS_IN_FragPos = world_pos.xyz;

	mat3 model_mat = mat3(model);

	PS_IN_Normal = normalize(model_mat * VS_IN_Normal);
    PS_IN_TexCoord = VS_IN_Texcoord;

	gl_Position = projection * view * world_pos;
}
