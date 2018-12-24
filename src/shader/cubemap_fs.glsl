
out vec4 PS_OUT_Color;

in vec3 PS_IN_TexCoord;

uniform samplerCube s_Skybox;
uniform vec3 white_point;
uniform float exposure;

void main()
{
	vec3 radiance = textureLod(s_Skybox, normalize(PS_IN_TexCoord), 0.0).rgb;

	radiance = pow(vec3(1,1,1) - exp(-radiance / white_point * exposure), vec3(1.0 / 2.2));

	PS_OUT_Color = vec4(radiance, 1.0);
}