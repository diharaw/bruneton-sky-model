#include "constants.glsl"
#include "uniforms.glsl"
#include "utility.glsl"
#include "transmittance_functions.glsl"
#include "irradiance_functions.glsl"

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x = LOCAL_SIZE, local_size_y = LOCAL_SIZE, local_size_z = 1) in;

// ------------------------------------------------------------------
// IMAGES -----------------------------------------------------------
// ------------------------------------------------------------------

layout (binding = 0, rgba32f) uniform image3D delta_irradiance;
layout (binding = 1, rgba32f) uniform image3D irradiance_read;
layout (binding = 2, rgba32f) uniform image3D irradiance_write;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

uniform vec4 blend;
uniform int scattering_order;

uniform sampler3D single_rayleigh_scattering;
uniform sampler3D single_mie_scattering;
uniform sampler3D multiple_scattering;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec2 coord = gl_GlobalInvocationID.xy;
    vec2 frag_coord = coord + vec2(0.5, 0.5);

    vec3 delta_irradiance = ComputeIndirectIrradianceTexture(single_rayleigh_scattering, single_mie_scattering, multiple_scattering, frag_coord, scattering_order);

    vec3 irradiance = RadianceToLuminance(delta_irradiance);

    imageStore(delta_irradiance, coord, vec4(delta_irradiance, 1.0));

    imageStore(irradiance_write, coord, vec4(irradiance, 1.0));

    if(blend[1] == 1)
        imageStore(irradiance_write, coord, imageLoad(irradiance_write, coord) + imageLoad(irradiance_read, coord));
}

// ------------------------------------------------------------------