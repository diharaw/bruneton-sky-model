// Using a fullscreen triangle for post-processing.
// https://www.saschawillems.de/?page_id=2122
// https://michaldrobot.com/2014/04/01/gcn-execution-patterns-in-full-screen-passes/

layout (std140) uniform u_GlobalUBO
{ 
    mat4 view;
    mat4 projection;
};

out vec3 PS_IN_TexCoord;

// void main()
// {
//     vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
// 	gl_Position = vec4(pos * 2.0f + -1.0f, 0.0f, 1.0f);

//     PS_IN_TexCoord = mat3(view) * vec3(pos, 1.0);
// }

void main(void)
{
     const vec3 vertices[4] = vec3[4](vec3(-1.0f, -1.0f, 1.0f),
                                      vec3( 1.0f, -1.0f, 1.0f),
                                      vec3(-1.0f,  1.0f, 1.0f),
                                      vec3( 1.0f,  1.0f, 1.0f));
    PS_IN_TexCoord = mat3(view) * vertices[gl_VertexID];
    gl_Position = vec4(vertices[gl_VertexID], 1.0f);
}
