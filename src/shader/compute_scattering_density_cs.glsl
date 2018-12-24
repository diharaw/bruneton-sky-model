#include <constants.glsl>
#include <uniforms.glsl>
#include <utility.glsl>
#include <transmittance_functions.glsl>
#include <scattering_functions.glsl>
#include <irradiance_functions.glsl>

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x = LOCAL_SIZE, local_size_y = LOCAL_SIZE, local_size_z = 1) in;

// ------------------------------------------------------------------
// IMAGES -----------------------------------------------------------
// ------------------------------------------------------------------

layout (binding = 0, rgba32f) uniform image3D delta_scattering_density;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

uniform vec4 blend;
uniform int layer;
uniform int scattering_order;

uniform sampler2D transmittance;
uniform sampler3D single_rayleigh_scattering;
uniform sampler3D single_mie_scattering;
uniform sampler3D multiple_scattering;
uniform sampler2D irradiance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    coord.z = layer;
    vec3 frag_coord = coord + vec3(0.5,0.5,0.5);

    vec3 scattering_density = ComputeScatteringDensityTexture(transmittance, single_rayleigh_scattering, single_mie_scattering, multiple_scattering, irradiance, frag_coord, scattering_order);

    imageStore(delta_scattering_density, coord, vec4(scattering_density, 1.0));
}

// ------------------------------------------------------------------