#include <uniforms.glsl>
#include <utility.glsl>
#include <transmittance_functions.glsl>
#include <scattering_functions.glsl>
#include <irradiance_functions.glsl>
#include <rendering_functions.glsl>

out vec4 PS_OUT_Color;

in vec2 PS_IN_TexCoord;
in vec3 PS_IN_FragPos;
in vec3 PS_IN_Normal;

// ------------------------------------------------------------------
// CONSTANTS --------------------------------------------------------
// ------------------------------------------------------------------

const vec3 kSphereCenter = vec3(0.0, 1.0, 0.0);
const float kSphereRadius = 1.0;
const vec3 kSphereAlbedo = vec3(0.8, 0.8, 0.8);
const vec3 kGroundAlbedo = vec3(0.04, 0.04, 0.04);

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

layout (std140) uniform u_GlobalUBO
{ 
    mat4 view;
    mat4 projection;
	mat4 inv_view;
	mat4 inv_projection;
	vec4 camera_pos;
};

uniform float exposure;
uniform vec3 color;
uniform vec3 white_point;
uniform vec3 earth_center;
uniform vec3 sun_direction;
uniform vec2 sun_size;

uniform sampler2D transmittance_texture;
uniform sampler2D irradiance_texture;
uniform sampler3D scattering_texture;
uniform sampler3D single_mie_scattering_texture;

// ------------------------------------------------------------------
// FUNCTIONS --------------------------------------------------------
// ------------------------------------------------------------------

/*
The functions to compute shadows and light shafts must be defined before we
can use them in the main shader function, so we define them first. Testing if
a point is in the shadow of the sphere S is equivalent to test if the
corresponding light ray intersects the sphere, which is very simple to do.
However, this is only valid for a punctual light source, which is not the case
of the Sun. In the following function we compute an approximate (and biased)
soft shadow by taking the angular size of the Sun into account:
*/
float GetSunVisibility(vec3 _point, vec3 sun_direction)
{
	vec3 p = _point - kSphereCenter;
	float p_dot_v = dot(p, sun_direction);
	float p_dot_p = dot(p, p);
	float ray_sphere_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
	float distance_to_intersection = -p_dot_v - sqrt(max(0.0, kSphereRadius * kSphereRadius - ray_sphere_center_squared_distance));
	if (distance_to_intersection > 0.0) 
	{
		// Compute the distance between the view ray and the sphere, and the
		// corresponding (tangent of the) subtended angle. Finally, use this to
		// compute an approximate sun visibility.
		float ray_sphere_distance = kSphereRadius - sqrt(ray_sphere_center_squared_distance);
		float ray_sphere_angular_distance = -ray_sphere_distance / p_dot_v;
		return smoothstep(1.0, 0.0, ray_sphere_angular_distance / sun_size.x);
	}
	return 1.0;
}

// ------------------------------------------------------------------

/*
The sphere also partially occludes the sky light, and we approximate this
effect with an ambient occlusion factor. The ambient occlusion factor due to a
sphere is given in <a href=
"http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf"
>Radiation View Factors</a> (Isidoro Martinez, 1995). In the simple case where
the sphere is fully visible, it is given by the following function:
*/
float GetSkyVisibility(vec3 _point) 
{
	vec3 p = _point - kSphereCenter;
	float p_dot_p = dot(p, p);
	return 1.0 + p.y / sqrt(p_dot_p) * kSphereRadius * kSphereRadius / p_dot_p;
}

// ------------------------------------------------------------------

/*
To compute light shafts we need the intersections of the view ray with the
shadow volume of the sphere S. Since the Sun is not a punctual light source this
shadow volume is not a cylinder but a cone (for the umbra, plus another cone for
the penumbra, but we ignore it here):
*/
void GetSphereShadowInOut(vec3 view_direction, vec3 sun_direction, out float d_in, out float d_out)
{
	vec3 camera = camera_pos.xyz;
	vec3 pos = camera - kSphereCenter;
	float pos_dot_sun = dot(pos, sun_direction);
	float view_dot_sun = dot(view_direction, sun_direction);
	float k = sun_size.x;
	float l = 1.0 + k * k;
	float a = 1.0 - l * view_dot_sun * view_dot_sun;
	float b = dot(pos, view_direction) - l * pos_dot_sun * view_dot_sun -
		k * kSphereRadius * view_dot_sun;
	float c = dot(pos, pos) - l * pos_dot_sun * pos_dot_sun -
		2.0 * k * kSphereRadius * pos_dot_sun - kSphereRadius * kSphereRadius;
	float discriminant = b * b - a * c;
	if (discriminant > 0.0) 
	{
		d_in = max(0.0, (-b - sqrt(discriminant)) / a);
		d_out = (-b + sqrt(discriminant)) / a;
		// The values of d for which delta is equal to 0 and kSphereRadius / k.
		float d_base = -pos_dot_sun / view_dot_sun;
		float d_apex = -(pos_dot_sun + kSphereRadius / k) / view_dot_sun;
		if (view_dot_sun > 0.0) 
		{
			d_in = max(d_in, d_apex);
			d_out = a > 0.0 ? min(d_out, d_base) : d_base;
		}
		else 
		{
			d_in = a > 0.0 ? max(d_in, d_base) : d_base;
			d_out = min(d_out, d_apex);
		}
	}
	else 
	{
		d_in = 0.0;
		d_out = 0.0;
	}
}

// ------------------------------------------------------------------

#ifdef RADIANCE_API_ENABLED

RadianceSpectrum GetSolarRadiance() 
{
	return solar_irradiance / (PI * sun_angular_radius * sun_angular_radius);
}

// ------------------------------------------------------------------

RadianceSpectrum GetSkyRadiance(
	Position camera, Direction view_ray, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	return GetSkyRadiance(transmittance_texture,
		scattering_texture, single_mie_scattering_texture,
		camera, view_ray, shadow_length, sun_direction, transmittance);
}

// ------------------------------------------------------------------

RadianceSpectrum GetSkyRadianceToPoint(
	Position camera, Position _point, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	return GetSkyRadianceToPoint(transmittance_texture,
		scattering_texture, single_mie_scattering_texture,
		camera, _point, shadow_length, sun_direction, transmittance);
}

// ------------------------------------------------------------------

IrradianceSpectrum GetSunAndSkyIrradiance(
	Position p, Direction normal, Direction sun_direction,
	out IrradianceSpectrum sky_irradiance) 
{
	return GetSunAndSkyIrradiance(transmittance_texture,
		irradiance_texture, p, normal, sun_direction, sky_irradiance);
}

// ------------------------------------------------------------------

#else

Luminance3 GetSolarRadiance()
{
	return solar_irradiance /
		(PI * sun_angular_radius * sun_angular_radius) *
		SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

// ------------------------------------------------------------------

Luminance3 GetSkyRadiance(
	Position camera, Direction view_ray, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	return GetSkyRadiance(transmittance_texture,
		scattering_texture, single_mie_scattering_texture,
		camera, view_ray, shadow_length, sun_direction, transmittance) *
		SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

// ------------------------------------------------------------------

Luminance3 GetSkyRadianceToPoint(
	Position camera, Position _point, Length shadow_length,
	Direction sun_direction, out DimensionlessSpectrum transmittance) 
{
	return GetSkyRadianceToPoint(transmittance_texture,
		scattering_texture, single_mie_scattering_texture,
		camera, _point, shadow_length, sun_direction, transmittance) *
		SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

// ------------------------------------------------------------------

Illuminance3 GetSunAndSkyIrradiance(
	Position p, Direction normal, Direction sun_direction,
	out IrradianceSpectrum sky_irradiance) 
{
	IrradianceSpectrum sun_irradiance = GetSunAndSkyIrradiance(
		transmittance_texture, irradiance_texture, p, normal,
		sun_direction, sky_irradiance);
	sky_irradiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
	return sun_irradiance * SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

#endif

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
	vec3 sky_irradiance;
	vec3 sun_irradiance = GetSunAndSkyIrradiance(PS_IN_FragPos - earth_center, PS_IN_Normal, sun_direction, sky_irradiance);

	vec3 sphere_radiance = color * (1.0 / PI) * (sun_irradiance + sky_irradiance);

	vec3 transmittance;
	vec3 in_scatter = GetSkyRadianceToPoint(camera_pos.xyz - earth_center, PS_IN_FragPos - earth_center, 0.0, sun_direction, transmittance);

	sphere_radiance = sphere_radiance * transmittance + in_scatter;

	sphere_radiance = pow(vec3(1,1,1) - exp(-sphere_radiance / white_point * exposure), vec3(1.0 / 2.2));

	PS_OUT_Color = vec4(sphere_radiance, 1.0);
}

// ------------------------------------------------------------------