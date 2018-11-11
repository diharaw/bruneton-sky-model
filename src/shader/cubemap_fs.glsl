
out vec4 PS_OUT_Color;

in vec3 PS_IN_TexCoord;

uniform samplerCube s_Skybox;

void main()
{
	PS_OUT_Color = vec4(texture(s_Skybox, normalize(PS_IN_TexCoord)).rgb, 1.0);
}
